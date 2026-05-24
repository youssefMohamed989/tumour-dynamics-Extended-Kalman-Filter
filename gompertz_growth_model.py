"""
gompertz_growth_model.py
========================
Generalized three-parameter Gompertz growth model with:
  - Curve fitting via nonlinear least squares
  - Bootstrap 95% confidence intervals on all parameters
  - Linear and exponential alternative model comparison
  - Prediction with uncertainty bands
  - Visualization

The Gompertz model is widely applicable wherever growth approaches
an asymptote in a sigmoid-with-skew manner:

    y(t) = A * exp(-B * exp(-C * t))

    A  = asymptotic (carrying-capacity) value
    B  = dimensionless displacement / shape parameter  (B > 0 shifts
         the inflection point rightward along the time axis)
    C  = intrinsic growth-rate constant  [units: 1 / time_unit]

Inflection point occurs at t* = ln(B) / C,  y* = A / e

Application examples
--------------------
Domain           y(t)                    t            Typical A
Oncology         Tumour volume (mm³)      Days         100–1000 mm³
Ecology          Population count         Years        carrying cap.
Microbiology     Colony / OD600           Hours        max OD
Economics        Cumulative revenue ($)   Quarters     market size
Epidemiology     Cumulative cases         Days         final plateau
Battery aging    Capacity retention (%)   Cycles       ~80 %
Product adoption Cumulative users         Months       TAM

Reference: Norton L., Cancer Res 1988;48:7067 (original oncology
           application); Tjørve & Tjørve, PLoS ONE 2017 (review).

Usage
-----
    from gompertz_growth_model import GompertzModel

    model = GompertzModel(time_unit="days", value_label="Tumour volume (mm³)")
    model.fit(t, y, y_sd=sd_array)
    model.bootstrap_ci(n_boot=1000)
    model.compare_alternatives()
    model.plot()
    model.summary()
    y_pred, lower, upper = model.predict(t_new)
"""

import numpy as np
from scipy.optimize import curve_fit, least_squares
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings


# ---------------------------------------------------------------------------
# Core model function (module-level, picklable for multiprocessing)
# ---------------------------------------------------------------------------

def _gompertz(t, A, B, C):
    """Three-parameter Gompertz:  y = A * exp(-B * exp(-C * t))"""
    return A * np.exp(-B * np.exp(-C * np.asarray(t, dtype=float)))


def _linear(t, m, b):
    return m * np.asarray(t, dtype=float) + b


def _exponential(t, y0, k):
    return y0 * np.exp(k * np.asarray(t, dtype=float))


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class GompertzModel:
    """
    Fit and interrogate a Gompertz growth model for any time-series dataset.

    Parameters
    ----------
    time_unit : str
        Label for the time axis (default "time").
    value_label : str
        Label for the y-axis (default "y").
    seed : int
        Random seed for reproducible bootstrap resampling.
    """

    def __init__(self, time_unit="time", value_label="y", seed=42):
        self.time_unit   = time_unit
        self.value_label = value_label
        self._rng        = np.random.default_rng(seed)

        # Populated after fit()
        self.t          = None
        self.y          = None
        self.y_sd       = None
        self.params     = None      # (A, B, C)
        self.pcov       = None
        self.residuals  = None
        self.r2         = None
        self.mape       = None

        # Populated after bootstrap_ci()
        self.ci         = None      # dict  {A:(lo,hi), B:(lo,hi), C:(lo,hi)}
        self.boot_params = None     # shape (n_boot, 3)

        # Populated after compare_alternatives()
        self.alt_results = {}

    # ------------------------------------------------------------------
    # Fitting
    # ------------------------------------------------------------------

    def fit(self, t, y, y_sd=None,
            A_init=None, B_init=5.0, C_init=0.1,
            bounds_A=(1e-3, np.inf),
            bounds_B=(1e-3, 100.0),
            bounds_C=(1e-5, 10.0)):
        """
        Fit the Gompertz model to data.

        Parameters
        ----------
        t, y : array-like
            Time points and observed values.
        y_sd : array-like or None
            Standard deviations for each observation (used as sigma
            weights in the fit and for error-bar plotting). Pass None
            to fit unweighted.
        A_init : float or None
            Initial guess for asymptote A.  Defaults to 1.5 * max(y).
        B_init, C_init : float
            Initial guesses for displacement B and growth rate C.
        bounds_A/B/C : 2-tuples
            (lower, upper) parameter bounds.
        """
        self.t    = np.asarray(t, dtype=float)
        self.y    = np.asarray(y, dtype=float)
        self.y_sd = np.asarray(y_sd, dtype=float) if y_sd is not None else None

        if A_init is None:
            A_init = max(self.y) * 1.5

        p0 = [A_init, B_init, C_init]
        lb = [bounds_A[0], bounds_B[0], bounds_C[0]]
        ub = [bounds_A[1], bounds_B[1], bounds_C[1]]
        sigma = self.y_sd if self.y_sd is not None else None

        try:
            popt, pcov = curve_fit(
                _gompertz, self.t, self.y,
                p0=p0, bounds=(lb, ub),
                sigma=sigma, absolute_sigma=True,
                maxfev=50_000
            )
        except RuntimeError as exc:
            raise RuntimeError(
                f"Gompertz fit did not converge. Try adjusting initial "
                f"guesses or bounds.\nOriginal error: {exc}"
            ) from exc

        self.params    = popt
        self.pcov      = pcov
        y_hat          = _gompertz(self.t, *popt)
        self.residuals = self.y - y_hat
        ss_res = np.sum(self.residuals ** 2)
        ss_tot = np.sum((self.y - self.y.mean()) ** 2)
        self.r2   = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        self.mape = float(np.mean(np.abs(self.residuals / self.y)) * 100)
        return self

    # ------------------------------------------------------------------
    # Bootstrap confidence intervals
    # ------------------------------------------------------------------

    def bootstrap_ci(self, n_boot=1000, alpha=0.05):
        """
        Residual-resampling bootstrap for 95% CIs on A, B, C.

        Parameters
        ----------
        n_boot : int
            Number of bootstrap resamples.
        alpha : float
            Significance level (default 0.05 → 95% CI).

        Returns
        -------
        ci : dict
            {'A': (lo, hi), 'B': (lo, hi), 'C': (lo, hi)}
        """
        if self.params is None:
            raise RuntimeError("Call fit() before bootstrap_ci().")

        y_fit  = _gompertz(self.t, *self.params)
        resid  = self.residuals
        boots  = np.zeros((n_boot, 3))

        for b in range(n_boot):
            resid_star = resid[self._rng.integers(0, len(resid), len(resid))]
            y_star     = y_fit + resid_star
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    popt_b, _ = curve_fit(
                        _gompertz, self.t, y_star,
                        p0=self.params,
                        bounds=([1e-3, 1e-3, 1e-5], [np.inf, 100.0, 10.0]),
                        maxfev=10_000
                    )
                boots[b] = popt_b
            except RuntimeError:
                boots[b] = self.params  # fall back to point estimate

        lo_pct = alpha / 2 * 100
        hi_pct = (1 - alpha / 2) * 100
        self.boot_params = boots
        self.ci = {
            "A": (np.percentile(boots[:, 0], lo_pct), np.percentile(boots[:, 0], hi_pct)),
            "B": (np.percentile(boots[:, 1], lo_pct), np.percentile(boots[:, 1], hi_pct)),
            "C": (np.percentile(boots[:, 2], lo_pct), np.percentile(boots[:, 2], hi_pct)),
        }
        return self.ci

    # ------------------------------------------------------------------
    # Alternative models
    # ------------------------------------------------------------------

    def compare_alternatives(self):
        """
        Fit linear and exponential alternatives and report R² and MAPE.
        Results stored in self.alt_results and printed.
        """
        if self.params is None:
            raise RuntimeError("Call fit() before compare_alternatives().")

        def _fit_alt(fn, p0, bounds=None):
            kw = dict(p0=p0, maxfev=50_000)
            if bounds:
                kw["bounds"] = bounds
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    popt, _ = curve_fit(fn, self.t, self.y, **kw)
            except RuntimeError:
                return None, np.nan, np.nan
            y_hat = fn(self.t, *popt)
            ss_res = np.sum((self.y - y_hat) ** 2)
            ss_tot = np.sum((self.y - self.y.mean()) ** 2)
            r2   = 1.0 - ss_res / ss_tot
            mape = float(np.mean(np.abs((self.y - y_hat) / self.y)) * 100)
            return popt, r2, mape

        p_lin, r2_lin, mape_lin = _fit_alt(_linear, [1.0, 0.0])
        p_exp, r2_exp, mape_exp = _fit_alt(
            _exponential, [self.y[0], 0.01],
            bounds=([0, 0], [np.inf, np.inf])
        )

        self.alt_results = {
            "linear":      {"params": p_lin, "r2": r2_lin, "mape": mape_lin},
            "exponential": {"params": p_exp, "r2": r2_exp, "mape": mape_exp},
        }
        return self.alt_results

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------

    def predict(self, t_new):
        """
        Predict y at new time points with bootstrap uncertainty bands.

        Parameters
        ----------
        t_new : array-like
            Time points at which to predict.

        Returns
        -------
        y_mean : ndarray   – point estimate from fitted parameters
        y_lower : ndarray  – 2.5th bootstrap percentile (or None if no CI)
        y_upper : ndarray  – 97.5th bootstrap percentile (or None if no CI)
        """
        if self.params is None:
            raise RuntimeError("Call fit() first.")
        t_new = np.asarray(t_new, dtype=float)
        y_mean = _gompertz(t_new, *self.params)

        if self.boot_params is not None:
            boot_preds = np.array([_gompertz(t_new, *p) for p in self.boot_params])
            y_lower = np.percentile(boot_preds, 2.5, axis=0)
            y_upper = np.percentile(boot_preds, 97.5, axis=0)
        else:
            y_lower = y_upper = None

        return y_mean, y_lower, y_upper

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self):
        """Print a formatted parameter table."""
        if self.params is None:
            print("Model not yet fitted.")
            return

        A, B, C = self.params
        t_star  = np.log(B) / C          # inflection time
        y_star  = A / np.e               # value at inflection

        print("=" * 55)
        print("  Gompertz Growth Model – Summary")
        print("=" * 55)
        print(f"  A  (asymptote)           = {A:.4g}")
        print(f"  B  (displacement)        = {B:.4f}")
        print(f"  C  (growth rate)         = {C:.6f}  [1/{self.time_unit}]")
        print(f"  Inflection time t*       = {t_star:.2f}  {self.time_unit}")
        print(f"  Inflection value y*      = {y_star:.4g}")
        print(f"  R²                       = {self.r2:.4f}")
        print(f"  MAPE                     = {self.mape:.2f}%")

        if self.ci:
            print("\n  Bootstrap 95% CIs:")
            for name, (lo, hi) in self.ci.items():
                print(f"    {name}  [{lo:.4g},  {hi:.4g}]")

        if self.alt_results:
            print("\n  Model comparison:")
            print(f"    {'Model':<14}  R²      MAPE")
            print(f"    {'Gompertz':<14}  {self.r2:.4f}  {self.mape:.2f}%")
            for key, res in self.alt_results.items():
                print(f"    {key.capitalize():<14}  {res['r2']:.4f}  {res['mape']:.2f}%")
        print("=" * 55)

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------

    def plot(self, t_min=None, t_max=None, figsize=(13, 9), save_path=None):
        """
        Four-panel diagnostic figure:
          (1) Fits + CI band, (2) Residual stem, (3) Bootstrap C dist,
          (4) Bootstrap A dist.
        """
        if self.params is None:
            raise RuntimeError("Call fit() before plot().")

        A, B, C = self.params
        t_lo = t_min if t_min is not None else self.t[0] * 0.9
        t_hi = t_max if t_max is not None else self.t[-1] * 1.1
        t_fine = np.linspace(t_lo, t_hi, 400)

        y_fine, y_lo_fine, y_hi_fine = self.predict(t_fine)

        fig = plt.figure(figsize=figsize)
        gs  = gridspec.GridSpec(2, 4, figure=fig, hspace=0.4, wspace=0.45)
        ax1 = fig.add_subplot(gs[0, :2])   # fits
        ax2 = fig.add_subplot(gs[0, 2:])   # residuals
        ax3 = fig.add_subplot(gs[1, :2])   # bootstrap C
        ax4 = fig.add_subplot(gs[1, 2:])   # bootstrap A

        # -- Panel 1: Fits --
        if y_lo_fine is not None:
            ax1.fill_between(t_fine, y_lo_fine, y_hi_fine,
                             alpha=0.25, color="#4A90D9", label="Gompertz 95% CI")
        ax1.plot(t_fine, y_fine, "-", color="#1A5FA8", lw=2.5,
                 label=f"Gompertz (R²={self.r2:.3f})")

        if self.alt_results:
            if self.alt_results["linear"]["params"] is not None:
                ax1.plot(t_fine, _linear(t_fine, *self.alt_results["linear"]["params"]),
                         "--", color="#C0541A", lw=1.8,
                         label=f"Linear (R²={self.alt_results['linear']['r2']:.3f})")
            if self.alt_results["exponential"]["params"] is not None:
                ax1.plot(t_fine, _exponential(t_fine, *self.alt_results["exponential"]["params"]),
                         "-.", color="#2E8B38", lw=1.8,
                         label=f"Exponential (R²={self.alt_results['exponential']['r2']:.3f})")

        if self.y_sd is not None:
            ax1.errorbar(self.t, self.y, yerr=1.96 * self.y_sd,
                         fmt="ko", ms=7, lw=1.5, capsize=4,
                         label="Observed (mean ± 1.96 SD)")
        else:
            ax1.plot(self.t, self.y, "ko", ms=7, label="Observed")

        # Mark inflection point
        t_star = np.log(B) / C
        y_star = A / np.e
        if t_lo <= t_star <= t_hi:
            ax1.axvline(t_star, color="purple", lw=1, ls=":", alpha=0.7)
            ax1.axhline(y_star, color="purple", lw=1, ls=":", alpha=0.7)
            ax1.scatter([t_star], [y_star], color="purple", zorder=5, s=60,
                        label=f"Inflection (t*={t_star:.1f})")

        ax1.set_xlabel(f"Time ({self.time_unit})", fontsize=10)
        ax1.set_ylabel(self.value_label, fontsize=10)
        ax1.set_title("Model Comparison", fontsize=11, fontweight="bold")
        ax1.legend(fontsize=7.5, loc="upper left")
        ax1.grid(True, alpha=0.3)

        # -- Panel 2: Residuals --
        ax2.stem(self.t, self.residuals, linefmt="C0-", markerfmt="C0o",
                 basefmt="k-")
        ax2.axhline(0, color="k", lw=1)
        ax2.set_xlabel(f"Time ({self.time_unit})", fontsize=10)
        ax2.set_ylabel(f"Residual ({self.value_label.split('(')[-1].rstrip(')')})"
                       if "(" in self.value_label else "Residual", fontsize=10)
        ax2.set_title(f"Gompertz Residuals\nMAPE = {self.mape:.2f}%",
                      fontsize=11, fontweight="bold")
        ax2.grid(True, alpha=0.3)

        # -- Panels 3 & 4: Bootstrap distributions --
        for ax, param_idx, name, color in [
            (ax3, 2, "C  (growth rate)", "#2A7FBF"),
            (ax4, 0, "A  (asymptote)",   "#BF6A2A"),
        ]:
            if self.boot_params is not None:
                vals = self.boot_params[:, param_idx]
                ax.hist(vals, bins=30, color=color, alpha=0.75, density=True,
                        edgecolor="none")
                point = self.params[param_idx]
                lo, hi = self.ci[name[0]]
                ax.axvline(point, color="red", lw=1.8,
                           label=f"Estimate = {point:.4g}")
                ax.axvline(lo,    color="k", lw=1, ls="--", label="95% CI")
                ax.axvline(hi,    color="k", lw=1, ls="--")
                ax.set_xlabel(f"Parameter {name}", fontsize=10)
                ax.set_ylabel("Density", fontsize=10)
                ax.set_title(f"Bootstrap Distribution\n{name[0]} = {point:.4g} "
                             f"[{lo:.4g}, {hi:.4g}]", fontsize=10, fontweight="bold")
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
            else:
                ax.text(0.5, 0.5, "Run bootstrap_ci() first",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=11, color="grey")

        fig.suptitle("Gompertz Growth Model – Diagnostic Summary",
                     fontsize=13, fontweight="bold", y=1.01)

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.tight_layout()
        plt.show()
        return fig


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------

def fit_gompertz(t, y, y_sd=None, n_boot=500, time_unit="time",
                 value_label="y", seed=42, plot=True, verbose=True):
    """
    One-call convenience wrapper.

    Returns
    -------
    model : GompertzModel   (fully fitted, with CI and alt-model comparison)

    Example
    -------
    >>> t  = [63, 77, 91, 112]               # days
    >>> y  = [32.1, 89.4, 178.2, 298.6]     # mm³
    >>> sd = [8.2,  18.3, 31.5,  47.2]
    >>> m  = fit_gompertz(t, y, y_sd=sd, time_unit="days",
    ...                   value_label="Tumour volume (mm³)")
    >>> print(m.predict([120, 130, 150]))    # extrapolate
    """
    m = GompertzModel(time_unit=time_unit, value_label=value_label, seed=seed)
    m.fit(t, y, y_sd=y_sd)
    m.bootstrap_ci(n_boot=n_boot)
    m.compare_alternatives()
    if verbose:
        m.summary()
    if plot:
        m.plot()
    return m


# ---------------------------------------------------------------------------
# Self-test / demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # --- Tumour volume (oncology) ---
    t_days = [63, 77, 91, 112]
    V_mean = [32.1, 89.4, 178.2, 298.6]
    V_sd   = [8.2,  18.3, 31.5,  47.2]

    print("\n=== Demo 1: Murine tumour volume ===")
    model = fit_gompertz(
        t_days, V_mean, y_sd=V_sd,
        n_boot=200, time_unit="days",
        value_label="Tumour volume (mm³)"
    )

    # --- Bacterial growth (OD600) ---
    print("\n=== Demo 2: Bacterial growth curve ===")
    t_hours = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]
    od600   = [0.02, 0.03, 0.05, 0.09, 0.18, 0.38, 0.65, 1.12,
               1.45, 1.61, 1.72, 1.75, 1.76]
    fit_gompertz(
        t_hours, od600, n_boot=200,
        time_unit="hours", value_label="OD600"
    )
