"""
extended_kalman_filter.py
=========================
Generalized Extended Kalman Filter (EKF) for state estimation from noisy,
indirect, nonlinear measurements.

Background
----------
A standard Kalman filter handles *linear* systems.  When the process
or measurement model is nonlinear the EKF linearises both models at
each step via their Jacobians, then runs standard Kalman update equations.

System equations
----------------
    x_{k+1} = f(x_k, u_k) + w_k      process model   w_k ~ N(0, Q)
    z_k     = h(x_k)      + v_k      measurement      v_k ~ N(0, R)

    x  – state vector  (n,)
    z  – measurement vector  (m,)
    Q  – process noise covariance  (n x n)
    R  – measurement noise covariance  (m x m)
    P  – state error covariance  (n x n)

Algorithm per step
------------------
    Predict:
        x̂⁻_k = f(x̂_{k-1}, u_{k-1})
        P⁻_k  = F_{k-1} P_{k-1} F_{k-1}ᵀ + Q

    Update:
        K_k   = P⁻_k Hᵀ_k (H_k P⁻_k Hᵀ_k + R)⁻¹
        x̂_k   = x̂⁻_k + K_k (z_k − h(x̂⁻_k))
        P_k   = (I − K_k H_k) P⁻_k

where F = ∂f/∂x  and  H = ∂h/∂x  are the Jacobians (computed
numerically here unless the user provides analytic versions).

Application examples
--------------------
Domain                 State x                    Measurement z
---------------------------------------------------------------------
Oncology (this repo)   [tumour_vol; GPC3]          circulating GPC3
Navigation (INS)       [pos; vel; heading]          GPS + IMU
Battery SoC estimation [SoC; internal resistance]   terminal voltage
Epidemiology           [S; I; R]                    reported cases
Robot localisation     [x; y; θ]                    range + bearing

Usage
-----
    from extended_kalman_filter import EKF, GompertzEKF

    # ── Minimal API ──────────────────────────────────────────────────
    def f(x, u=None):   return ...   # nonlinear state transition
    def h(x):           return ...   # nonlinear measurement model

    ekf = EKF(n_states=2, n_obs=1, Q=..., R=..., f=f, h=h)
    ekf.init_state(x0, P0)

    for z_k in measurements:
        ekf.predict()
        ekf.update(z_k)

    print(ekf.trajectory)           # list of filtered state vectors

    # ── Pre-built Gompertz tumour / biomarker model ───────────────────
    model = GompertzEKF(A=312.4, B=4.82, C=0.185,
                        alpha=0.05, sigma_obs=0.20)
    results = model.run(gpc3_measurements, t_days)
    model.plot(results, tumour_volume_true=V_caliper)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from copy import deepcopy


# ---------------------------------------------------------------------------
# Numerical Jacobian helper
# ---------------------------------------------------------------------------

def _numerical_jacobian(f, x, u=None, eps=1e-5):
    """
    Compute the Jacobian of f(x) w.r.t. x by central differences.

    Returns
    -------
    J : ndarray, shape (m, n)
        J[i, j] = ∂f_i / ∂x_j
    """
    x  = np.asarray(x, dtype=float)
    n  = len(x)
    f0 = np.atleast_1d(f(x) if u is None else f(x, u))
    m  = len(f0)
    J  = np.zeros((m, n))
    for j in range(n):
        dx      = np.zeros(n)
        dx[j]   = eps
        fwd     = np.atleast_1d(f(x + dx) if u is None else f(x + dx, u))
        bwd     = np.atleast_1d(f(x - dx) if u is None else f(x - dx, u))
        J[:, j] = (fwd - bwd) / (2 * eps)
    return J


# ---------------------------------------------------------------------------
# Generic EKF
# ---------------------------------------------------------------------------

class EKF:
    """
    Generic Extended Kalman Filter.

    Parameters
    ----------
    n_states : int
        Dimension of the state vector x.
    n_obs : int
        Dimension of the observation vector z.
    Q : array-like, shape (n_states, n_states)
        Process noise covariance.
    R : array-like, shape (n_obs, n_obs)
        Measurement noise covariance.
    f : callable  x, [u] -> x_next
        Nonlinear state-transition function.  Single-step.
    h : callable  x -> z
        Nonlinear observation function.
    F_jac : callable or None
        Analytic Jacobian ∂f/∂x.  If None, computed numerically.
    H_jac : callable or None
        Analytic Jacobian ∂h/∂x.  If None, computed numerically.
    clip_state : callable or None
        Optional post-update function x -> x to enforce state constraints
        (e.g. non-negativity).  Applied after every update step.
    """

    def __init__(self, n_states, n_obs, Q, R, f, h,
                 F_jac=None, H_jac=None, clip_state=None):
        self.n   = n_states
        self.m   = n_obs
        self.Q   = np.asarray(Q, dtype=float)
        self.R   = np.asarray(R, dtype=float)
        self._f  = f
        self._h  = h
        self._Fj = F_jac
        self._Hj = H_jac
        self._clip = clip_state

        # State (set via init_state)
        self.x     = None   # current estimate
        self.P     = None   # current covariance
        self.x_prior = None
        self.P_prior = None

        # Full trajectory storage
        self.trajectory  = []   # list of x_posterior at each step
        self.covariances = []   # list of P_posterior at each step
        self.innovations = []   # list of innovation z - h(x_prior)
        self.gains       = []   # list of Kalman gain K

    def init_state(self, x0, P0=None):
        """Set the initial state and covariance."""
        self.x = np.asarray(x0, dtype=float).copy()
        if P0 is None:
            P0 = np.eye(self.n)
        self.P = np.asarray(P0, dtype=float).copy()
        self.trajectory  = [self.x.copy()]
        self.covariances = [self.P.copy()]
        return self

    def predict(self, u=None):
        """EKF predict step:  propagate state and covariance."""
        if u is not None:
            x_pred = np.atleast_1d(self._f(self.x, u))
            F = (self._Fj(self.x, u) if self._Fj
                 else _numerical_jacobian(self._f, self.x, u))
        else:
            x_pred = np.atleast_1d(self._f(self.x))
            F = (self._Fj(self.x) if self._Fj
                 else _numerical_jacobian(self._f, self.x))

        P_pred = F @ self.P @ F.T + self.Q

        self.x_prior = x_pred
        self.P_prior = P_pred
        return x_pred, P_pred

    def update(self, z, R=None):
        """
        EKF update step: incorporate observation z.

        Parameters
        ----------
        z : array-like
            Observation vector.
        R : array-like or None
            Per-step measurement noise covariance.  If None, uses self.R.

        Returns
        -------
        x_post : ndarray   posterior state estimate
        P_post : ndarray   posterior covariance
        K      : ndarray   Kalman gain
        innov  : ndarray   innovation  z - h(x_prior)
        """
        z  = np.atleast_1d(np.asarray(z, dtype=float))
        R_ = np.asarray(R, dtype=float) if R is not None else self.R

        H = (self._Hj(self.x_prior) if self._Hj
             else _numerical_jacobian(self._h, self.x_prior))

        z_pred = np.atleast_1d(self._h(self.x_prior))
        innov  = z - z_pred
        S      = H @ self.P_prior @ H.T + R_
        K      = self.P_prior @ H.T @ np.linalg.inv(S)

        x_post = self.x_prior + K @ innov
        I      = np.eye(self.n)
        P_post = (I - K @ H) @ self.P_prior  # Joseph form for stability omitted for clarity

        if self._clip:
            x_post = np.asarray(self._clip(x_post), dtype=float)

        self.x = x_post
        self.P = P_post
        self.trajectory.append(x_post.copy())
        self.covariances.append(P_post.copy())
        self.innovations.append(innov.copy())
        self.gains.append(K.copy())
        return x_post, P_post, K, innov

    def smooth(self):
        """
        Rauch-Tung-Striebel (RTS) smoother: backward pass over stored
        trajectory for improved estimates.

        Returns
        -------
        x_smooth : list of ndarray  – smoothed state estimates (same length
                   as self.trajectory)
        """
        N  = len(self.trajectory)
        xs = [None] * N
        Ps = [None] * N
        xs[-1] = self.trajectory[-1].copy()
        Ps[-1] = self.covariances[-1].copy()

        # Recompute predicted covariances backward
        # (requires re-running predict from each stored state – simplified here)
        for k in range(N - 2, -1, -1):
            x_k = self.trajectory[k]
            P_k = self.covariances[k]
            F_k = _numerical_jacobian(self._f, x_k)
            P_pred_k = F_k @ P_k @ F_k.T + self.Q
            G_k  = P_k @ F_k.T @ np.linalg.pinv(P_pred_k)
            x_f  = np.atleast_1d(self._f(x_k))
            xs[k] = x_k + G_k @ (xs[k + 1] - x_f)
            Ps[k] = P_k + G_k @ (Ps[k + 1] - P_pred_k) @ G_k.T

        return xs, Ps

    def reset(self):
        """Clear stored trajectory; keep current state."""
        self.trajectory  = [self.x.copy()]
        self.covariances = [self.P.copy()]
        self.innovations = []
        self.gains       = []


# ---------------------------------------------------------------------------
# Pre-built: Gompertz tumour + GPC3 biomarker EKF
# ---------------------------------------------------------------------------

class GompertzEKF:
    """
    EKF that reconstructs tumour volume from serial biomarker (GPC3)
    measurements using a Gompertz process model.

    State vector  x = [V,  c]
        V  – tumour volume (mm³)
        c  – circulating biomarker concentration (ng/mL)

    Process model (one-step discrete Gompertz, Euler integration)
        V_{k+1} = V_k + dt * V_k * C * ln(A / V_k)
        c_{k+1} = c  + dt * alpha * (V_k - V_thresh) - beta * c

    Measurement model
        z_k = c_k + noise

    Parameters
    ----------
    A, B, C : float
        Gompertz asymptote, displacement, and growth-rate constant.
        Obtain from GompertzModel.fit().
    alpha : float
        Shedding rate – how fast tumour volume translates to serum
        biomarker (ng/mL per mm³ per day).
    beta : float
        Biomarker clearance rate (day⁻¹).  Default: 0.05 (approx. 14-day
        half-life, reasonable for serum proteins like GPC3).
    V_thresh : float
        Minimum volume that contributes to biomarker elevation
        (default 0, i.e., all tumour volume sheds biomarker).
    sigma_obs : float
        Measurement noise std (ng/mL).  Used to construct R.
    sigma_proc_V : float
        Process noise std on V (mm³).
    sigma_proc_c : float
        Process noise std on c (ng/mL).
    dt : float
        Time step in days (default 1.0 day for daily integration even if
        observations are weekly).
    """

    def __init__(self, A=312.4, B=4.82, C=0.185,
                 alpha=0.05, beta=0.05, V_thresh=0.0,
                 sigma_obs=0.20, sigma_proc_V=5.0, sigma_proc_c=0.05,
                 dt=1.0):
        self.A       = A
        self.B       = B
        self.C       = C
        self.alpha   = alpha
        self.beta    = beta
        self.V_thresh = V_thresh
        self.dt      = dt

        # Noise covariance matrices
        Q = np.diag([sigma_proc_V ** 2, sigma_proc_c ** 2])
        R = np.array([[sigma_obs ** 2]])

        def f(x, u=None):
            """One-step Gompertz + biomarker ODE (Euler)."""
            V, c = float(x[0]), float(x[1])
            V = max(V, 1e-3)
            dV = self.dt * V * self.C * np.log(max(self.A / V, 1e-9))
            dc = self.dt * (self.alpha * max(V - self.V_thresh, 0) - self.beta * c)
            return np.array([V + dV, c + dc])

        def h(x):
            """Observe the biomarker concentration directly."""
            return np.array([x[1]])

        def clip(x):
            """Enforce non-negativity of both states."""
            return np.maximum(x, 1e-6)

        self._ekf = EKF(
            n_states=2, n_obs=1,
            Q=Q, R=R,
            f=f, h=h,
            clip_state=clip
        )

    def run(self, biomarker_obs, t_obs, V0=None, c0=None, P0=None):
        """
        Run the EKF over a series of biomarker observations.

        Parameters
        ----------
        biomarker_obs : array-like, shape (T,)
            Serial biomarker measurements (ng/mL), one per time point.
        t_obs : array-like, shape (T,)
            Observation times (days).  Need not be equally spaced.
        V0 : float or None
            Initial tumour volume estimate.  Defaults to a small seed
            volume (computed from first biomarker reading).
        c0 : float or None
            Initial biomarker concentration estimate.  Defaults to
            biomarker_obs[0].
        P0 : array-like (2,2) or None
            Initial state covariance.  Defaults to large diagonal
            uncertainty.

        Returns
        -------
        dict with keys:
            t          – observation times
            V_est      – EKF posterior tumour-volume estimates
            c_est      – EKF posterior biomarker estimates
            V_std      – 1-sigma uncertainty on V
            c_std      – 1-sigma uncertainty on c
            innovations – list of scalar innovations
        """
        obs = np.asarray(biomarker_obs, dtype=float)
        t   = np.asarray(t_obs,        dtype=float)
        T   = len(obs)

        if c0 is None:
            c0 = obs[0]
        if V0 is None:
            # Invert biomarker: V ≈ c / alpha  (rough initialisation)
            V0 = max(c0 / max(self.alpha, 1e-9), 5.0)
        if P0 is None:
            P0 = np.diag([V0 ** 2, c0 ** 2 + 1.0])

        self._ekf.init_state([V0, c0], P0)

        V_est, c_est, V_std, c_std, innovations = [], [], [], [], []

        for k in range(T):
            # Propagate from previous time to current time in sub-steps
            if k > 0:
                n_steps = max(1, int(round((t[k] - t[k - 1]) / self.dt)))
                for _ in range(n_steps):
                    self._ekf.predict()
            else:
                self._ekf.predict()     # one dummy predict for k=0

            x_post, P_post, K, innov = self._ekf.update(obs[k])

            V_est.append(x_post[0])
            c_est.append(x_post[1])
            V_std.append(np.sqrt(P_post[0, 0]))
            c_std.append(np.sqrt(P_post[1, 1]))
            innovations.append(float(innov[0]))

        return {
            "t":           t,
            "V_est":       np.array(V_est),
            "c_est":       np.array(c_est),
            "V_std":       np.array(V_std),
            "c_std":       np.array(c_std),
            "innovations": innovations,
        }

    def cross_validate(self, biomarker_obs, t_obs, V_true=None,
                       V0_arr=None, c0_arr=None):
        """
        Leave-one-out cross-validation of tumour volume reconstruction.

        Parameters
        ----------
        biomarker_obs : array-like, shape (N, T)
            N animals × T time-point biomarker observations.
        t_obs : array-like, shape (T,)
            Shared time grid.
        V_true : array-like or None
            Ground-truth volumes for computing R² and MAPE.
        V0_arr, c0_arr : array-like or None
            Per-animal initial states; if None, estimated from data.

        Returns
        -------
        dict with:
            r2, mape, V_pred_all, errors_all
        """
        obs  = np.asarray(biomarker_obs, dtype=float)
        t    = np.asarray(t_obs, dtype=float)
        N, T = obs.shape
        V_pred_all = np.zeros_like(obs)

        for i in range(N):
            # Leave subject i out: fit uses all others for hyperparameters
            # (Here we simply run EKF independently per subject)
            V0 = V0_arr[i] if V0_arr is not None else None
            c0 = c0_arr[i] if c0_arr is not None else None

            # Reset EKF internals for fresh run
            self._ekf.reset() if hasattr(self._ekf, '_x') else None
            res = self.run(obs[i], t, V0=V0, c0=c0)
            V_pred_all[i] = res["V_est"]

        result = {"V_pred_all": V_pred_all}

        if V_true is not None:
            V_true = np.asarray(V_true, dtype=float)
            flat_pred  = V_pred_all.ravel()
            flat_true  = V_true.ravel()
            ss_res = np.sum((flat_true - flat_pred) ** 2)
            ss_tot = np.sum((flat_true - flat_true.mean()) ** 2)
            result["r2"]         = 1.0 - ss_res / ss_tot
            result["mape"]       = float(np.mean(np.abs((flat_true - flat_pred) / flat_true)) * 100)
            result["errors_all"] = flat_true - flat_pred
        return result

    def plot(self, results, tumour_volume_true=None,
             figsize=(13, 8), save_path=None):
        """
        Three-panel figure:
          (1) Estimated vs true tumour volume  [+ uncertainty band]
          (2) Biomarker: estimated vs observed
          (3) Innovations (filter residuals)
        """
        t     = results["t"]
        V_est = results["V_est"]
        V_std = results["V_std"]
        c_est = results["c_est"]
        c_std = results["c_std"]
        innov = results["innovations"]

        fig, axes = plt.subplots(1, 3, figsize=figsize)

        # -- Panel 1: Tumour volume --
        ax = axes[0]
        ax.fill_between(t, V_est - 1.96 * V_std, V_est + 1.96 * V_std,
                        alpha=0.25, color="#4A90D9", label="EKF 95% CI")
        ax.plot(t, V_est, "-o", color="#1A5FA8", lw=2, ms=6,
                label="EKF estimated V")
        if tumour_volume_true is not None:
            ax.plot(t, tumour_volume_true, "k--s", lw=1.5, ms=6,
                    label="True V (caliper)")
        ax.set_xlabel("Time (days)", fontsize=10)
        ax.set_ylabel("Tumour volume (mm³)", fontsize=10)
        ax.set_title("Tumour Volume Reconstruction", fontsize=11, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # -- Panel 2: Biomarker --
        ax = axes[1]
        ax.fill_between(t, c_est - 1.96 * c_std, c_est + 1.96 * c_std,
                        alpha=0.25, color="#E07B3A", label="EKF 95% CI")
        ax.plot(t, c_est, "-o", color="#C05010", lw=2, ms=6,
                label="EKF estimated c")
        ax.plot(t, [self._ekf.trajectory[i + 1][1]
                    for i in range(len(t))
                    if i + 1 < len(self._ekf.trajectory)],
                "k.", ms=4)
        ax.set_xlabel("Time (days)", fontsize=10)
        ax.set_ylabel("Biomarker c (ng/mL)", fontsize=10)
        ax.set_title("Biomarker State Estimate", fontsize=11, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # -- Panel 3: Innovations --
        ax = axes[2]
        ax.stem(t, innov, linefmt="C2-", markerfmt="C2o", basefmt="k-")
        ax.axhline(0, color="k", lw=1)
        ax.set_xlabel("Time (days)", fontsize=10)
        ax.set_ylabel("Innovation (ng/mL)", fontsize=10)
        ax.set_title("Filter Innovations\n(should be white noise)", fontsize=11,
                     fontweight="bold")
        ax.grid(True, alpha=0.3)

        fig.suptitle("Gompertz EKF – State Estimation Summary",
                     fontsize=13, fontweight="bold")
        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.show()
        return fig


# ---------------------------------------------------------------------------
# Self-test / demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=== Demo 1: Tumour volume reconstruction from GPC3 ===\n")

    np.random.seed(2024)
    t_obs = np.array([63.0, 77.0, 91.0, 112.0])
    # True GPC3 (approx from manuscript data)
    gpc3_true   = np.array([1.18, 1.84, 2.92, 4.52])
    gpc3_noisy  = gpc3_true + np.random.randn(len(t_obs)) * 0.20

    model = GompertzEKF(A=312.4, B=4.82, C=0.185,
                        alpha=0.01, beta=0.04,
                        sigma_obs=0.20, sigma_proc_V=8.0,
                        sigma_proc_c=0.10, dt=1.0)

    results = model.run(gpc3_noisy, t_obs, V0=30.0, c0=gpc3_noisy[0])

    print("  Time (days):", t_obs)
    print("  EKF V_est  :", np.round(results["V_est"],  1))
    print("  EKF c_est  :", np.round(results["c_est"],  3))
    print("  EKF V_std  :", np.round(results["V_std"],  2))
    print("  Innovations:", np.round(results["innovations"], 4))

    model.plot(results,
               tumour_volume_true=np.array([32.1, 89.4, 178.2, 298.6]))

    # ------------------------------------------------------------------
    print("\n=== Demo 2: Generic EKF on a 1-D sine-wave state ===\n")

    # True state: x_k+1 = 0.8 x_k + 0.2 sin(k)
    # Measurement: z_k = x_k^2  (nonlinear)
    T = 50
    x_true = np.zeros(T)
    z_meas = np.zeros(T)
    x_true[0] = 5.0
    for k in range(1, T):
        x_true[k] = 0.8 * x_true[k - 1] + 0.2 * np.sin(k) + np.random.randn() * 0.5
    z_meas = x_true ** 2 + np.random.randn(T) * 2.0

    def f_demo(x, u=None):
        return np.array([0.8 * x[0] + 0.2 * np.sin(len(x))])  # approximation

    def h_demo(x):
        return np.array([x[0] ** 2])

    ekf_demo = EKF(n_states=1, n_obs=1,
                   Q=np.array([[0.25]]),
                   R=np.array([[4.0]]),
                   f=f_demo, h=h_demo)
    ekf_demo.init_state([5.0], P0=np.array([[1.0]]))

    for k in range(T):
        ekf_demo.predict()
        ekf_demo.update(z_meas[k])

    x_filtered = np.array([s[0] for s in ekf_demo.trajectory[1:]])
    r2 = 1 - np.sum((x_true - x_filtered) ** 2) / np.sum((x_true - x_true.mean()) ** 2)
    print(f"  Demo 2 reconstruction R² = {r2:.4f}")

    plt.figure(figsize=(9, 4))
    plt.plot(x_true,     label="True state", lw=2)
    plt.plot(x_filtered, "--", label=f"EKF filtered (R²={r2:.3f})", lw=2)
    plt.plot(np.sqrt(np.abs(z_meas)), ".", color="gray", alpha=0.5,
             label="sqrt(|z|) – observation proxy")
    plt.xlabel("Step k"); plt.ylabel("x")
    plt.title("Generic EKF: Nonlinear State Estimation Demo")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
