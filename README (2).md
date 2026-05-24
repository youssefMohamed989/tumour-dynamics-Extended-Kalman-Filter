# SiNP-Chitosan Nanoplasmonic Lab-on-a-Chip for Whole-Blood GPC3 Detection

A self-powered, 3D-printed microfluidic biosensor that performs passive plasma separation and quantitative glypican-3 (GPC3) detection directly from unprocessed finger-prick whole blood, without centrifugation, pipettes, cold-chain reagents, or any laboratory instrument.

---

## Overview

This repository contains all custom code accompanying the manuscript:

**"A SiNP-Chitosan Nanocomposite Converts Hydrophobic 3D-Printed Microchannels into Whole-Blood-Operable Nanoplasmonic Biosensors"**

Youssef M. Hassan, Mohamed S. Attia, Hala El-Tantawi, Ibrahim Rabie Ali, Dalia M. El-Husseini, Mona N. Abou-Omar

The platform integrates four functional modules in a single 45 x 18 x 6 mm SLA-printed cartridge:

1. Geometry-programmed capillary transport driven by a SiNP-chitosan nanocomposite surface coating
2. Micropillar-assisted passive plasma separation (35-42% recovery in 30 s)
3. Anti-GPC3-functionalised AuNP nanoplasmonic biosensing (LOD = 0.24 ng/mL, LSPR channel)
4. Dual-mode readout: LSPR ratiometry (A580/A520) for laboratory validation and instrument-free RGB colourimetry (TCS3200/Arduino) for point-of-care use

Clinical validation was performed in a DEN-induced murine HCC model across weeks 9-16 (n = 10 HCC, n = 10 control). All performance metrics are preclinical; prospective human validation is planned.

---

## Repository Structure

```
.
|-- arduino/
|   `-- rgb_sensor_readout/
|       `-- rgb_sensor_readout.ino     # TCS3200 RGB sensor + Arduino Nano firmware
|
|-- matlab/
|   |-- surface_dynamics/
|   |   `-- surface_dynamics.m         # Surface property modelling (Ra, theta, Laplace pressure, fibrinogen)
|   |-- cfd_validation/
|   |   `-- carreau_yasuda_flow.m      # Non-Newtonian capillary flow model and CY vs Newtonian comparison
|   |-- biosensor_calibration/
|   |   `-- lspr_rgb_calibration.m    # LSPR ratiometric and RGB calibration curve fitting
|   |-- tumour_dynamics/
|   |   |-- gompertz_fit.m             # Gompertz tumour growth model fitting
|   |   `-- ekf_tumour_reconstruction.m  # Extended Kalman Filter for tumour volume reconstruction from GPC3
|   `-- statistics/
|       `-- biomarker_statistics.m     # ROC analysis, Bland-Altman, Passing-Bablok, PCA
|
`-- README.md
```

---

## Hardware

### Chip Fabrication

- Printer: Formlabs Form 3, Clear V4 Resin
- Layer height: 25 um
- Post-cure: 405 nm UV, 60 degrees C, 30 min (FormCure)
- Chip dimensions: 45 x 18 x 6 mm
- Micropillar diameter: 50-100 um; inter-pillar gap: 150 um; height ~0.45 mm
- Design software: Autodesk Fusion 360 (STL files available upon request from corresponding author)

### Surface Coating (SiNP-Chitosan Nanocomposite)

Sequential deposition:

1. Stober-synthesised SiNPs (~142 nm): spin-coat at 2000 rpm for 60 s, anneal at 80 degrees C for 20 min
2. TPP-crosslinked chitosan (2% w/v, 50-190 kDa, >=85% deacetylation, pH 5.5): spin-coat at 1500 rpm for 45 s, dry at 37 degrees C for 2 h

Resulting surface properties:
- Contact angle: 74.3 degrees (bare resin) to 25.3 degrees (nanocomposite)
- Laplace pressure: 86 Pa to 451 Pa (5.2-fold amplification)
- Fibrinogen adsorption: 3.8 ug/cm2 to 0.31 ug/cm2 (12-fold reduction)

### AuNP Biosensor Assembly

- AuNP synthesis: Frens-Turkevich method; citrate-stabilised; TEM core diameter 30.2 +/- 2.8 nm; zeta = -38.2 mV (pre-conjugation); PDI = 0.09
- Antibody conjugation: EDC/Sulfo-NHS carbodiimide coupling; anti-GPC3 (Abcam ab56789, clone 1G12); 10 ug/mL in PBS pH 7.4; 30 min room temperature
- Blocking: 1% BSA overnight; stored in PBS + 0.05% sodium azide at 4 degrees C
- Post-conjugation characterisation: hydrodynamic diameter 45.2 +/- 3.1 nm (PDI 0.11); zeta = -26.4 +/- 1.8 mV
- Stable for 45 days at 37 degrees C: A580/A520 drift <3%, DLS diameter change <4 nm

### Readout Electronics

- RGB channel (point-of-care): TCS3200 colour sensor + Arduino Nano
- LSPR channel (laboratory validation): Jenway UV-1800 UV-Vis spectrophotometer; A580/A520 ratiometry
- Arduino IDE version: 2.3.2

---

## Software

### Requirements

#### MATLAB

- MATLAB R2023b
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- Signal Processing Toolbox (for EKF)

#### Arduino

- Arduino IDE 2.3.2
- Board: Arduino Nano (ATmega328P)
- No external libraries required beyond standard Arduino core

#### Statistical Analysis (R)

All statistical analyses reported in the manuscript were performed in R v4.3.1. The MATLAB scripts reproduce the key figures and models; for full statistical output including Benjamini-Hochberg FDR-corrected q-values, refer to Supplementary Table S-E18.

---

## Module Descriptions

### arduino/rgb_sensor_readout/rgb_sensor_readout.ino

Firmware for the TCS3200 colour sensor module integrated into the chip detection chamber.

Reads raw red, green, and blue channel photodiode counts from the TCS3200, applies a simple calibration offset, and outputs a formatted serial string at 9600 baud for data logging. The red channel intensity is the primary GPC3 quantification signal (calibration: y = -4.54x + 116.68; R2 = 0.994; LOD = 0.31 ng/mL).

Wiring:
- S0 -> D4, S1 -> D5 (frequency scaling, set to 20%)
- S2 -> D6, S3 -> D7 (colour filter selection)
- OUT -> D8 (frequency output)
- OE -> GND (output enable, active low)

Serial output format: `R:<value> G:<value> B:<value> GPC3_est:<value_ng_per_mL>`

### matlab/surface_dynamics/surface_dynamics.m

Reproduces Figure 4 of the manuscript. Computes and plots the progressive surface property changes across the three surface conditions (bare SLA resin, SiNP-only, SiNP-chitosan nanocomposite) using experimentally measured inputs:

- Surface roughness Ra from profilometry
- Contact angle theta from goniometry
- Laplace capillary pressure from theta and channel geometry
- Fibrinogen adsorption from QCM-D
- Coating thickness from ellipsometry and cross-sectional SEM

Inputs are defined as constants at the top of the script matching the values reported in Table 1 and Supporting Information Tables S-E7 through S-E9. Model fidelity R2 = 0.982 +/- 0.009 (n = 18) vs experimental profilometry.

### matlab/cfd_validation/carreau_yasuda_flow.m

Implements the Carreau-Yasuda non-Newtonian capillary flow model and compares it against the Newtonian assumption and experimental flow-front position data (n = 10 chips). Reproduces Figure 3.

Carreau-Yasuda parameters used:
- eta_0 = 0.056 Pa.s (zero-shear viscosity)
- eta_inf = 0.0035 Pa.s (infinite-shear viscosity)
- lambda = 3.313 s (relaxation time)
- a = 2.0, n = 0.3568 (power-law index and transition parameter)

Key results reproduced:
- CY model: RMSE = 7.3%, R2 = 0.994 vs experiment
- Newtonian model: RMSE = 23.1%, R2 = 0.871 vs experiment
- Mean effective viscosity eta_eff = 8.9 mPa.s in the operating shear-rate window (20-50 s^-1)
- Wall shear stress: mean 0.94 Pa, peak 1.83 Pa (below haemolysis threshold of ~4 Pa)
- Theoretical flow rate: 6.4 uL/min (experimental: 6.1-6.8 uL/min)

### matlab/biosensor_calibration/lspr_rgb_calibration.m

Fits and plots calibration curves for both detection channels across 0.25-5 ng/mL recombinant GPC3:

- LSPR channel: A580/A520 ratiometric signal; multi-day replicated (n = 21 blank replicates across 3 assay days); R2 = 0.997; slope = 0.183 AU per ng/mL; LOD = 0.24 ng/mL; LOQ = 0.80 ng/mL
- RGB channel: TCS3200 red-channel intensity; R2 = 0.994; LOD = 0.31 ng/mL

LOD is calculated as mean blank + 3 sigma (n = 21 blank replicates). LOQ is calculated as mean blank + 10 sigma. Both channels are plotted with 95% confidence bands. Inter-channel Pearson correlation is computed from the paired calibration data (r = 0.994 in the manuscript).

### matlab/tumour_dynamics/gompertz_fit.m

Fits a three-parameter Gompertz growth model to the longitudinal tumour volume data from the DEN-HCC murine cohort (weeks 9-16; n = 10 HCC animals):

V(t) = A * exp(-B * exp(-C * t))

Fitted parameters: A = 312.4 mm3, B = 4.82, C = 0.185 day^-1 (manuscript values; bootstrap 95% CIs in Supplementary Table S8). Goodness of fit: R2 = 0.988, MAPE = 6.3%. Compares against linear (R2 = 0.82) and exponential (R2 = 0.91) alternatives.

Note: This model was fitted to the 10-animal HCC cohort. Independent prospective validation in a separate cohort is required before this constitutes a validated predictive tool.

### matlab/tumour_dynamics/ekf_tumour_reconstruction.m

Implements an Extended Kalman Filter (EKF) that reconstructs tumour volume trajectories from serial weekly GPC3 measurements alone, using the Gompertz state-transition equations as the process model.

State vector: [tumour_volume; GPC3_concentration]
Measurement: circulating GPC3 (on-chip LSPR)

Performance (leave-one-out cross-validated on the 10-animal fitting cohort): R2 = 0.960 (95% CI: 0.931-0.981). This is an in-sample retrospective proof-of-concept result. Independent validation in a separate animal cohort is required before clinical application.

### matlab/statistics/biomarker_statistics.m

Reproduces the statistical analyses reported in Sections 3.6-3.7:

- ROC analysis and AUC computation with DeLong method 95% CIs
- Bland-Altman agreement analysis between on-chip and ELISA GPC3 measurements
- Passing-Bablok regression (on-chip vs ELISA)
- Intraclass correlation coefficient (ICC, two-way mixed, absolute agreement)
- Principal component analysis (PCA) of the multi-biomarker panel (GPC3, AFP, ALT, AST, GGT)
- One-way ANOVA with Tukey HSD post-hoc correction
- Benjamini-Hochberg FDR correction for multi-timepoint comparisons

Input data format: CSV files with columns [animal_id, group, week, GPC3_chip, GPC3_elisa, AFP, ALT, AST, GGT, tumour_volume]. Data are available from the corresponding author upon reasonable request.

---

## Calibration

### LSPR Channel

Calibration standards: 0, 0.25, 0.5, 1, 2, 5 ng/mL recombinant murine GPC3 (R&D Systems cat. 4609-GP) prepared in normal mouse serum matrix. Run in triplicate on each assay day.

```
A580/A520 = 0.183 * [GPC3 ng/mL] + 0.312
LOD = 0.24 ng/mL
LOQ = 0.80 ng/mL
Linear range: 0.25 - 5 ng/mL
```

### RGB Channel

```
Red_intensity = -4.54 * [GPC3 ng/mL] + 116.68
LOD = 0.31 ng/mL
Linear range: 0.25 - 5 ng/mL
```

Both equations are from multi-day replicated calibration (R2 = 0.997 and 0.994 respectively). For single-run calibration the R2 will be lower due to run-to-run colloidal variance inherent in multivalent AuNP aggregation assays (representative single-run values shown in Figure 5D-E of the manuscript).

---

## Assay Protocol (Brief)

1. Load 5-8 uL whole blood (finger-prick or retro-orbital; heparin, EDTA, or citrate anticoagulated) at the chip inlet.
2. Passive capillary flow drives plasma separation across the SiNP-chitosan micropillar array within 30 s. No user intervention required.
3. Separated plasma contacts anti-GPC3-AuNP chitosan hydrogel beads in the detection zone.
4. Incubate 10-12 min at room temperature.
5. Read signal:
   - Point-of-care: red channel intensity from TCS3200/Arduino; apply RGB calibration equation.
   - Laboratory: measure A580 and A520 on benchtop UV-Vis; compute A580/A520 ratio; apply LSPR calibration equation.
6. Total sample-to-answer time: under 15 min.

Anticoagulant compatibility confirmed: heparin, EDTA, citrate (recovery >=95.8% for all). Haematocrit range validated: Hct 0.30-0.50 (no significant effect on GPC3 accuracy, ANOVA p = 0.22). Samples with spectrophotometric evidence of haemolysis (A414 > 0.05) should be excluded.

---

## Performance Summary

| Parameter | LSPR Channel | RGB Channel |
|---|---|---|
| LOD | 0.24 ng/mL | 0.31 ng/mL |
| LOQ | 0.80 ng/mL | ~1.0 ng/mL |
| Linear range | 0.25-5 ng/mL | 0.25-5 ng/mL |
| Multi-day R2 | 0.997 | 0.994 |
| Intra-day CV | <6% | <6% |
| Inter-day CV | <9% | <9% |
| Spike recovery | 97.3 +/- 2.1% | 97.3 +/- 2.1% |
| ELISA agreement (ICC) | 0.983 (95% CI: 0.962-0.994) | - |
| Inter-channel Pearson r | 0.994 | 0.994 |
| Conjugate stability | 45 days at 37 degrees C | 45 days at 37 degrees C |
| Assay time | <15 min | <15 min |

Murine DEN-HCC model (n = 10 HCC, n = 10 control; weeks 9-16):
- GPC3 pooled AUC: 1.00 (95% CI: 0.993-1.000)
- Tumour volume correlation: Pearson r = 0.981 (R2 = 0.98)

These values reflect ceiling-level performance in a tightly controlled, single-strain, single-sex, small-sample (n = 10 per group) murine cohort. They are not projections of human diagnostic performance. AUC values from controlled animal models carry a well-documented risk of optimistic bias due to the elimination of inter-individual biological heterogeneity. Prospective human validation across diverse HCC aetiologies is required before clinical deployment.

---

## Planned Validation Roadmap

- Phase 1 (6-12 months): Independent replication in a second BALB/c DEN-HCC cohort (n >=5 per group) at a separate laboratory.
- Phase 2 (12-18 months): Bridging feasibility study in archived de-identified human serum (n = 30 per group: confirmed HCC, liver cirrhosis without HCC, healthy controls).
- Phase 3 (24-36 months): Prospective multi-centre clinical validation across viral hepatitis B/C, NAFLD, and alcohol-related liver disease aetiologies.

Concurrent platform extensions: smartphone-based RGB image analysis; multiplexed GPC3-AFP-DCP panel on a single chip.

---

## Animal Ethics

All animal procedures were approved by the Institutional Animal Care and Use Committee (IACUC) of Ain Shams University, Faculty of Science (Approval Code: ASU-SCI/ZOOL/2024/2/4) and conducted in accordance with ARRIVE 2.0 guidelines. Twenty adult male BALB/c mice (6-8 weeks, 20-25 g) were obtained from the Theodor Bilharz Research Institute, Giza, Egypt.

---

## Data and Design File Availability

Raw experimental data, STL fabrication files, and additional fabrication parameters are available from the corresponding author upon reasonable request:

Youssef M. Hassan  
Department of Zoology, Faculty of Science, Ain Shams University, Abbassia 11566, Cairo, Egypt  
ORCID: 0009-0005-3615-4137  
Email: yousefmohamed_p@sci.asu.edu.eg

---

## Citation

If you use this code or platform design in your work, please cite the associated manuscript (citation details to be updated upon publication).

---

## Funding

This research received no external funding.

## Declaration of Competing Interests

The authors declare no known competing financial interests or personal relationships that could have appeared to influence the work reported.
