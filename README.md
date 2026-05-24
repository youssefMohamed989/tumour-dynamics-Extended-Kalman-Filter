# tumour-dynamics-Extended-Kalman-Filter

Gompertz tumour growth modelling and Extended Kalman Filter (EKF) state reconstruction of hepatocellular carcinoma burden from serial nanoplasmonic GPC3 biosensor measurements. Validated in a DEN-induced murine HCC model (n = 10/group, weeks 9–16). Includes non-Newtonian CFD, surface dynamics modelling, dual-mode LSPR/RGB calibration, and full biomarker statistics.

---

## Contents

```
.
├── gompertz_growth_model.py       # Generalised Gompertz fitting class (Python)
├── extended_kalman_filter.py      # Generic EKF + GompertzEKF pre-built model (Python)
├── gompertz_fit.m                 # Gompertz fitting + bootstrap CIs (MATLAB)
├── biomarker_statistics.m         # ROC, Bland-Altman, PCA, ANOVA, BH-FDR (MATLAB)
├── carreau_yasuda_flow.m          # Non-Newtonian capillary CFD (MATLAB)
├── lspr_rgb_calibration.m         # Dual-mode biosensor calibration curves (MATLAB)
├── surface_dynamics.m             # Surface property modelling (MATLAB)
└── rgb_sensor_readout.ino         # TCS3200 + Arduino Nano firmware (C++)
```

---

## Quick Start

### Python (recommended)

```bash
pip install numpy scipy matplotlib
```

```python
from gompertz_growth_model import fit_gompertz
from extended_kalman_filter import GompertzEKF

# 1. Fit Gompertz model to longitudinal volume data
t_days = [63, 77, 91, 112]                      # days post-DEN
V_mean = [32.1, 89.4, 178.2, 298.6]             # mm³ group mean
V_sd   = [8.2,  18.3, 31.5,  47.2]              # mm³ SD  (n = 10)

model = fit_gompertz(t_days, V_mean, y_sd=V_sd,
                     n_boot=1000,
                     time_unit="days",
                     value_label="Tumour volume (mm³)")
# → prints A, B, C with 95% CIs; saves 4-panel diagnostic figure

# 2. Run EKF to reconstruct volume from serial GPC3 alone
gpc3_obs = [1.18, 1.84, 2.92, 4.52]             # ng/mL (weekly chip readings)

ekf = GompertzEKF(A=model.params[0],
                  B=model.params[1],
                  C=model.params[2],
                  alpha=0.01,    # shedding rate  ng·mL⁻¹·mm⁻³·day⁻¹
                  beta=0.04,     # clearance rate day⁻¹
                  sigma_obs=0.20)

results = ekf.run(gpc3_obs, t_days)
print(results["V_est"])   # reconstructed tumour volumes
ekf.plot(results, tumour_volume_true=V_mean)
```

### MATLAB

```matlab
% Requires: Optimization Toolbox, Statistics & ML Toolbox, Signal Processing Toolbox
run('gompertz_fit.m')             % Gompertz fit + bootstrap
run('biomarker_statistics.m')     % ROC / Bland-Altman / PCA
run('carreau_yasuda_flow.m')      % CFD validation
run('lspr_rgb_calibration.m')     % Calibration curves
run('surface_dynamics.m')         % Surface property plots
```

---

## Python Modules

### `gompertz_growth_model.py`

Generalised three-parameter Gompertz model — no domain-specific hard-coding.

**Model:**
```
y(t) = A · exp(−B · exp(−C · t))

A  asymptotic (carrying-capacity) value
B  dimensionless displacement  →  inflection at  t* = ln(B)/C,  y* = A/e
C  intrinsic growth-rate constant  [1/time_unit]
```

**`GompertzModel` class — key methods:**

| Method | Description |
|---|---|
| `.fit(t, y, y_sd=None)` | Nonlinear least-squares fit via `scipy.optimize.curve_fit`; weighted if `y_sd` provided |
| `.bootstrap_ci(n_boot=1000)` | Residual-resampling bootstrap; returns `{'A':(lo,hi), 'B':…, 'C':…}` |
| `.compare_alternatives()` | Fits linear and exponential alternatives; reports R² and MAPE |
| `.predict(t_new)` | Returns `(y_mean, y_lower, y_upper)` with bootstrap uncertainty bands |
| `.summary()` | Prints parameter table, inflection point, goodness-of-fit |
| `.plot()` | Four-panel figure: fits + CI, residuals, bootstrap distributions for A and C |

**Convenience wrapper:**
```python
model = fit_gompertz(t, y, y_sd=sd, n_boot=500,
                     time_unit="hours", value_label="OD600")
y_pred, lo, hi = model.predict(np.linspace(0, 30, 200))
```

**Application domains:**

| Domain | `y(t)` | `t` | Typical `A` |
|---|---|---|---|
| Oncology | Tumour volume (mm³) | Days | 100–1000 mm³ |
| Microbiology | OD600 / colony count | Hours | Max OD |
| Ecology | Population size | Years | Carrying capacity |
| Epidemiology | Cumulative cases | Days | Final plateau |
| Battery aging | Capacity retention (%) | Cycles | ~80% |
| Product adoption | Cumulative users | Months | TAM |

---

### `extended_kalman_filter.py`

Two-layer design: a **generic EKF** that accepts arbitrary `f` and `h` callables, and a **pre-built `GompertzEKF`** wired to the Gompertz process model.

#### Generic `EKF` class

```python
from extended_kalman_filter import EKF

ekf = EKF(
    n_states = 2,
    n_obs    = 1,
    Q = np.diag([25.0, 0.01]),   # process noise covariance
    R = np.array([[0.04]]),       # measurement noise covariance
    f = my_process_fn,            # x_{k+1} = f(x_k)   [nonlinear]
    h = my_obs_fn,                # z_k      = h(x_k)   [nonlinear]
    # F_jac / H_jac: optional analytic Jacobians; numeric otherwise
    clip_state = lambda x: np.maximum(x, 0)   # enforce constraints
)

ekf.init_state(x0=[30.0, 1.0], P0=np.diag([900.0, 1.0]))

for z in measurements:
    ekf.predict()
    x_post, P_post, K, innov = ekf.update(z)

xs_smooth, Ps_smooth = ekf.smooth()   # RTS backward pass
```

Jacobians are computed by **central-difference numerical differentiation** unless analytic versions are supplied — making the class applicable to any nonlinear system without symbolic math.

#### State-space equations

```
Predict:
    x̂⁻_k  =  f(x̂_{k-1})
    P⁻_k   =  F_{k-1} P_{k-1} Fᵀ_{k-1} + Q

Update:
    K_k    =  P⁻_k Hᵀ_k (H_k P⁻_k Hᵀ_k + R)⁻¹
    x̂_k    =  x̂⁻_k + K_k (z_k − h(x̂⁻_k))
    P_k    =  (I − K_k H_k) P⁻_k

where  F = ∂f/∂x,  H = ∂h/∂x  (Jacobians)
```

#### Pre-built `GompertzEKF`

State vector `x = [V, c]`:

```
V_{k+1} = V_k + dt · V_k · C · ln(A / V_k)          (Gompertz ODE, Euler)
c_{k+1} = c_k + dt · (α · max(V_k − V_thresh, 0) − β · c_k)

z_k = c_k + noise                                     (measure biomarker only)
```

Key parameters:

| Parameter | Meaning | Default |
|---|---|---|
| `A, B, C` | Gompertz shape params (from `GompertzModel.fit`) | 312.4, 4.82, 0.185 |
| `alpha` | Biomarker shedding rate (ng·mL⁻¹·mm⁻³·day⁻¹) | 0.05 |
| `beta` | Biomarker serum clearance rate (day⁻¹) | 0.05 |
| `sigma_obs` | Measurement noise std (ng/mL) | 0.20 |
| `sigma_proc_V` | Process noise on V (mm³) | 5.0 |
| `sigma_proc_c` | Process noise on c (ng/mL) | 0.05 |
| `dt` | Integration sub-step (days) | 1.0 |

**`GompertzEKF` methods:**

| Method | Returns |
|---|---|
| `.run(biomarker_obs, t_obs)` | `dict` with `V_est, c_est, V_std, c_std, innovations` |
| `.cross_validate(obs_matrix, t_obs, V_true)` | `r2, mape, V_pred_all, errors_all` |
| `.plot(results, tumour_volume_true=None)` | Three-panel figure: V trajectory, biomarker, innovations |

**Porting to a different domain:**

Only three things change — provide new `f` and `h` callables and retune `Q`, `R`:

```python
# Example: battery state-of-charge estimation
def f_battery(x, u=None):
    SoC, R_int = x
    # discharge model ...
    return np.array([SoC_next, R_int_next])

def h_battery(x):
    SoC, R_int = x
    V_terminal = OCV(SoC) - I_load * R_int
    return np.array([V_terminal])

ekf = EKF(n_states=2, n_obs=1,
          Q=np.diag([1e-4, 1e-6]),
          R=np.array([[1e-3]]),
          f=f_battery, h=h_battery)
```

Other direct substitutions: INS navigation `[pos; vel; heading]` + GPS/IMU, SIR epidemiology `[S; I; R]` + reported case counts, robot localisation `[x; y; θ]` + lidar range-bearing.

---

## MATLAB Modules

### `gompertz_fit.m`

Fits `V(t) = A·exp(−B·exp(−C·t))` to longitudinal tumour volume data using `lsqcurvefit` (Optimization Toolbox). Compares against linear and exponential alternatives. Bootstrap 95% CIs computed by residual resampling (n = 1000, `rng(2024)`). Produces a three-panel figure: model comparison, residuals, bootstrap distribution of C.

Manuscript values: `A = 312.4 mm³, B = 4.82, C = 0.185 day⁻¹, R² = 0.988, MAPE = 6.3%`.

### `biomarker_statistics.m`

Full statistical pipeline for the murine DEN-HCC cohort (Sections 3.6–3.7):

| Analysis | Function / Method |
|---|---|
| ROC + AUC | `perfcurve` + bootstrap DeLong CI |
| Method agreement | Bland-Altman bias, ±1.96 SD LoA, proportional bias test |
| Method comparison | Passing-Bablok regression (all pairwise slopes → median) |
| Reliability | ICC(2,1) two-way mixed, absolute agreement, F-dist 95% CI |
| Dimensionality | PCA on standardised 5-biomarker panel |
| Group comparison | One-way ANOVA + Tukey HSD (`anova1`, `multcompare`) |
| Multiple testing | Benjamini-Hochberg FDR step-up correction (20 tests) |

Input: CSV with columns `[animal_id, group, week, GPC3_chip, GPC3_elisa, AFP, ALT, AST, GGT, tumour_volume]`. Falls back to manuscript summary statistics if no CSV is present.

### `carreau_yasuda_flow.m`

Implements the Carreau-Yasuda viscosity model for whole blood:

```
η(γ̇) = η_∞ + (η_0 − η_∞) · (1 + (λγ̇)^a)^((n−1)/a)

η_0 = 0.056 Pa·s,  η_∞ = 0.0035 Pa·s,  λ = 3.313 s,  a = 2.0,  n = 0.3568
```

Integrates the modified Lucas-Washburn equation numerically (Euler, `dt = 0.01 s`) to compute flow-front position vs. time in a 45 mm microchannel. Compares CY vs. Newtonian assumption against experimental data (n = 10 chips). Produces Figure 3: flow-front, parity plot, Bland-Altman, and viscosity profile.

Key outputs: `CY R² = 0.994, RMSE = 7.3%`; `Newtonian R² = 0.871, RMSE = 23.1%`; mean η_eff = 8.9 mPa·s; wall shear stress 0.94 Pa (peak 1.83 Pa, below haemolysis threshold).

### `lspr_rgb_calibration.m`

Fits linear calibration curves for both readout channels over 0.25–5 ng/mL recombinant murine GPC3 (n = 3 replicates × 3 assay days):

```
LSPR:  A₅₈₀/A₅₂₀ = 0.183 · [GPC3] + 0.312    R² = 0.997   LOD = 0.24 ng/mL
RGB:   I_red       = −4.54 · [GPC3] + 116.68   R² = 0.994   LOD = 0.31 ng/mL

LOD  =  μ_blank + 3σ_blank   (n = 21 blank replicates, 3 assay days)
LOQ  =  μ_blank + 10σ_blank
```

Plots calibration curves with 95% prediction bands (t-distribution, n − 2 df). Computes inter-channel Pearson r from paired calibration fits.

### `surface_dynamics.m`

Computes and visualises surface property changes across three conditions (bare SLA resin → SiNP-only → SiNP-chitosan) from experimentally measured inputs:

- Ra (profilometry), θ (goniometry), ΔP_Laplace = −2γ cos θ / r_h, fibrinogen adsorption (QCM-D), coating thickness (ellipsometry)
- Produces Figure 4: bar charts, simulated 2D topography maps, and normalised performance heatmap

Model fidelity R² = 0.982 ± 0.009 (n = 18) vs. profilometry.

---

## Arduino Firmware

### `rgb_sensor_readout.ino`

TCS3200 colour sensor readout for instrument-free GPC3 quantification.

**Wiring:**

```
TCS3200  →  Arduino Nano
S0       →  D4   (frequency scaling: S0=HIGH, S1=LOW → 20%)
S1       →  D5
S2       →  D6   (filter select)
S3       →  D7
OUT      →  D8   (pulse-frequency output)
OE       →  GND  (active low, always enabled)
```

**Serial output** (9600 baud):
```
R:<counts> G:<counts> B:<counts> GPC3_est:<ng/mL>
```

**Calibration equation applied on-device:**
```cpp
const float SLOPE     = -4.54f;   // counts per ng/mL
const float INTERCEPT = 116.68f;  // counts
// [GPC3] = (I_red − INTERCEPT) / SLOPE
```

Readings are averaged over `N_AVERAGE = 5` acquisitions per reported value. Values below `LOD = 0.31 ng/mL` are reported as `<LOD`. Update `SLOPE` and `INTERCEPT` from a daily single-point recalibration against a known standard.

---

## Dependencies

| Environment | Packages / Toolboxes |
|---|---|
| Python ≥ 3.9 | `numpy`, `scipy`, `matplotlib` |
| MATLAB R2023b | Statistics and ML Toolbox, Optimization Toolbox, Signal Processing Toolbox |
| Arduino | IDE 2.3.2; board: Arduino Nano (ATmega328P); no external libraries |

Install Python dependencies:
```bash
pip install numpy scipy matplotlib
```

---

## Performance

| Metric | Value |
|---|---|
| LSPR LOD | 0.24 ng/mL |
| RGB LOD | 0.31 ng/mL |
| LSPR multi-day R² | 0.997 |
| RGB multi-day R² | 0.994 |
| ELISA agreement (ICC) | 0.983 (95% CI: 0.962–0.994) |
| Assay time | < 15 min |
| Gompertz fit R² | 0.988 (MAPE 6.3%) |
| EKF reconstruction R² | 0.960 (95% CI: 0.931–0.981, LOO-CV) |
| CY CFD R² | 0.994 (RMSE 7.3%) |

> **Scope note.** All figures reflect a tightly controlled, single-strain, single-sex murine cohort (n = 10/group). The AUC = 1.00 and EKF R² = 0.960 are in-sample retrospective results and are not projections of human clinical performance. Independent validation in a separate cohort is required before any clinical application.

---

## Validation Roadmap

| Phase | Timeline | Milestone |
|---|---|---|
| 1 | 6–12 months | Independent replication, second BALB/c DEN-HCC cohort (n ≥ 5/group) |
| 2 | 12–18 months | Bridging study, archived human serum (n = 30/group: HCC / cirrhosis / healthy) |
| 3 | 24–36 months | Prospective multi-centre trial (HBV, HCV, NAFLD, ALD aetiologies) |

Concurrent extensions: smartphone RGB image analysis; multiplexed GPC3-AFP-DCP panel on a single chip.

---

## Citation

Citation details to be updated upon manuscript publication. If you use this code, please cite the associated manuscript (Hassan et al., in preparation).

---

## Ethics

Animal procedures approved by the IACUC of Ain Shams University, Faculty of Science (Approval Code: ASU-SCI/ZOOL/2024/2/4), conducted in accordance with ARRIVE 2.0 guidelines.

## Funding

No external funding received.

## Contact

**Youssef M. Hassan** — Department of Zoology, Faculty of Science, Ain Shams University, Abbassia 11566, Cairo, Egypt  
ORCID: [0009-0005-3615-4137](https://orcid.org/0009-0005-3615-4137)  
Email: yousefmohamed_p@sci.asu.edu.eg
