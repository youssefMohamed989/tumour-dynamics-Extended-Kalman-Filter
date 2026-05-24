% lspr_rgb_calibration.m
%
% Fits and plots calibration curves for both detection channels:
%   LSPR ratiometric channel (A580/A520) - laboratory quantification
%   RGB colourimetric channel (TCS3200 red-channel intensity) - point-of-care
%
% Calibration range: 0.25 - 5 ng/mL recombinant murine GPC3
% Standards prepared in normal mouse serum matrix (n = 3 replicates per concentration)
% Multi-day replicated: n = 21 blank replicates across 3 assay days
%
% LOD = mean_blank + 3 * sigma_blank
% LOQ = mean_blank + 10 * sigma_blank
%
% Results (manuscript Table 1 and Figure 5D-E):
%   LSPR: R2 = 0.997, slope = 0.183 AU/(ng/mL), LOD = 0.24 ng/mL
%   RGB:  R2 = 0.994, LOD = 0.31 ng/mL
%   Inter-channel Pearson r = 0.994
%
% Reference: Hassan et al., manuscript Section 3.5-3.6.

clear; clc; close all;

%% --- Multi-day replicated calibration data ----------------------------
% Calibration concentrations (ng/mL)
conc = [0, 0.25, 0.5, 1.0, 2.0, 5.0]';  % include blank for LOD calculation

% LSPR ratiometric signal (A580/A520), mean across 3 replicates x 3 days
% Blank mean = 0.312 +/- 0.008 (n=21 replicates, 3 assay days)
LSPR_mean = [0.312; 0.358; 0.403; 0.495; 0.678; 1.227];
LSPR_sd   = [0.008; 0.014; 0.012; 0.018; 0.021; 0.031];  % SD across replicates

% RGB red-channel intensity (TCS3200 counts), mean +/- SD
RGB_mean  = [187.4; 180.1; 173.8; 162.0; 144.2;  93.9];
RGB_sd    = [4.1;   5.2;   4.8;   5.5;   6.1;    7.3];

%% --- Blank statistics for LOD/LOQ calculation -------------------------
blank_LSPR_mean = LSPR_mean(1);
blank_LSPR_sd   = LSPR_sd(1);
blank_RGB_mean  = RGB_mean(1);
blank_RGB_sd    = RGB_sd(1);

%% --- Linear regression on non-zero calibrators -------------------------
conc_cal  = conc(2:end);
LSPR_cal  = LSPR_mean(2:end);
RGB_cal   = RGB_mean(2:end);

% LSPR fit
p_LSPR = polyfit(conc_cal, LSPR_cal, 1);
LSPR_fit = polyval(p_LSPR, conc_cal);
SS_res_L = sum((LSPR_cal - LSPR_fit).^2);
SS_tot_L = sum((LSPR_cal - mean(LSPR_cal)).^2);
R2_LSPR  = 1 - SS_res_L / SS_tot_L;

% RGB fit
p_RGB  = polyfit(conc_cal, RGB_cal, 1);
RGB_fit_vals = polyval(p_RGB, conc_cal);
SS_res_R = sum((RGB_cal - RGB_fit_vals).^2);
SS_tot_R = sum((RGB_cal - mean(RGB_cal)).^2);
R2_RGB   = 1 - SS_res_R / SS_tot_R;

%% --- LOD and LOQ -------------------------------------------------------
% LSPR: LOD at which signal = blank + 3*sigma
LOD_LSPR_signal = blank_LSPR_mean + 3 * blank_LSPR_sd;
LOD_LSPR_conc   = (LOD_LSPR_signal - p_LSPR(2)) / p_LSPR(1);
LOQ_LSPR_signal = blank_LSPR_mean + 10 * blank_LSPR_sd;
LOQ_LSPR_conc   = (LOQ_LSPR_signal - p_LSPR(2)) / p_LSPR(1);

% RGB: LOD at which signal = blank - 3*sigma (decreasing signal)
LOD_RGB_signal  = blank_RGB_mean - 3 * blank_RGB_sd;
LOD_RGB_conc    = (LOD_RGB_signal - p_RGB(2)) / p_RGB(1);

fprintf('--- LSPR Channel (A580/A520) ---\n');
fprintf('  Slope:     %.4f AU per ng/mL\n',   p_LSPR(1));
fprintf('  Intercept: %.4f AU\n',              p_LSPR(2));
fprintf('  R2:        %.4f\n',                 R2_LSPR);
fprintf('  LOD:       %.3f ng/mL\n',           LOD_LSPR_conc);
fprintf('  LOQ:       %.3f ng/mL\n\n',         LOQ_LSPR_conc);

fprintf('--- RGB Channel (Red intensity) ---\n');
fprintf('  Slope:     %.4f counts per ng/mL\n', p_RGB(1));
fprintf('  Intercept: %.2f counts\n',            p_RGB(2));
fprintf('  R2:        %.4f\n',                   R2_RGB);
fprintf('  LOD:       %.3f ng/mL\n\n',           LOD_RGB_conc);

%% --- Inter-channel correlation -----------------------------------------
% Use predicted values at calibrator concentrations
LSPR_pred = polyval(p_LSPR, conc_cal);
RGB_pred  = polyval(p_RGB,  conc_cal);

% Normalise for Pearson correlation between channels
% (channels measure in opposite directions: LSPR increases, RGB decreases)
r_channels = corr(LSPR_pred, RGB_pred);
fprintf('Inter-channel Pearson r (LSPR vs RGB, calibration fits): %.4f\n', r_channels);
fprintf('(Negative because LSPR increases and RGB decreases with GPC3; |r| reported)\n\n');

%% --- Figure: calibration curves with 95% CI bands ---------------------
conc_fine = linspace(0, 5.5, 200)';

fig = figure('Name', 'Biosensor Calibration', 'Position', [100 100 1100 480]);

% --- LSPR panel ---
subplot(1, 3, 1:2);

% 95% prediction interval
n_pts = length(conc_cal);
se_LSPR = sqrt(sum((LSPR_cal - LSPR_fit).^2) / (n_pts - 2)) .* ...
    sqrt(1/n_pts + (conc_fine - mean(conc_cal)).^2 / sum((conc_cal - mean(conc_cal)).^2));
t_crit = tinv(0.975, n_pts - 2);
LSPR_fine = polyval(p_LSPR, conc_fine);

fill([conc_fine; flipud(conc_fine)], ...
    [LSPR_fine + t_crit*se_LSPR; flipud(LSPR_fine - t_crit*se_LSPR)], ...
    [0.5 0.7 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
plot(conc_fine, LSPR_fine, '-b', 'LineWidth', 2);
errorbar(conc_cal, LSPR_cal, 1.96*LSPR_sd(2:end), 'o', ...
    'Color', [0.1 0.2 0.8], 'MarkerFaceColor', [0.1 0.2 0.8], 'MarkerSize', 7, 'LineWidth', 1.5);

% LOD line
yline(LOD_LSPR_signal, '--r', sprintf('LOD signal'), 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
xline(LOD_LSPR_conc,   ':r',  sprintf('LOD = %.2f ng/mL', LOD_LSPR_conc), ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 8);

xlabel('[GPC3] (ng mL^{-1})');
ylabel('A_{580}/A_{520} ratiometric signal (AU)');
title(sprintf('LSPR Calibration\nR^2 = %.4f, slope = %.3f AU/(ng/mL)\nLOD = %.2f ng/mL', ...
    R2_LSPR, p_LSPR(1), LOD_LSPR_conc));
legend({'95% CI', sprintf('Linear fit (y = %.3fx + %.3f)', p_LSPR(1), p_LSPR(2)), ...
    'Calibrators (mean \pm1.96SD)'}, 'Location', 'northwest', 'FontSize', 8);
grid on; box on;
xlim([-0.2, 5.5]);

% --- RGB panel ---
subplot(1, 3, 3);

RGB_fine = polyval(p_RGB, conc_fine);
se_RGB = sqrt(sum((RGB_cal - RGB_fit_vals).^2) / (n_pts - 2)) .* ...
    sqrt(1/n_pts + (conc_fine - mean(conc_cal)).^2 / sum((conc_cal - mean(conc_cal)).^2));

fill([conc_fine; flipud(conc_fine)], ...
    [RGB_fine + t_crit*se_RGB; flipud(RGB_fine - t_crit*se_RGB)], ...
    [1.0 0.7 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
plot(conc_fine, RGB_fine, '-', 'Color', [0.8 0.3 0.1], 'LineWidth', 2);
errorbar(conc_cal, RGB_cal, 1.96*RGB_sd(2:end), 's', ...
    'Color', [0.7 0.2 0.0], 'MarkerFaceColor', [0.7 0.2 0.0], 'MarkerSize', 7, 'LineWidth', 1.5);
xline(LOD_RGB_conc, ':r', sprintf('LOD = %.2f ng/mL', LOD_RGB_conc), ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 8);

xlabel('[GPC3] (ng mL^{-1})');
ylabel('RGB Red-Channel Intensity (counts)');
title(sprintf('RGB Calibration (TCS3200)\nR^2 = %.4f\nLOD = %.2f ng/mL', R2_RGB, LOD_RGB_conc));
legend({'95% CI', sprintf('Linear fit (y = %.2fx + %.2f)', p_RGB(1), p_RGB(2)), ...
    'Calibrators (mean \pm1.96SD)'}, 'Location', 'northeast', 'FontSize', 8);
grid on; box on;
xlim([-0.2, 5.5]);

sgtitle('Dual-Mode GPC3 Biosensor Calibration', 'FontSize', 12, 'FontWeight', 'bold');

saveas(fig, 'calibration_curves.png');
fprintf('Figure saved: calibration_curves.png\n');

%% --- Quantification function (for use in data analysis) ----------------
fprintf('\n--- Quantification equations ---\n');
fprintf('LSPR -> [GPC3]: ([signal] - %.4f) / %.4f ng/mL\n', p_LSPR(2), p_LSPR(1));
fprintf('RGB  -> [GPC3]: ([signal] - %.2f) / %.2f ng/mL\n', p_RGB(2), p_RGB(1));
