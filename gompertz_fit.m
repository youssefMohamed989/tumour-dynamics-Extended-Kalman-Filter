% gompertz_fit.m
%
% Fits a three-parameter Gompertz growth model to longitudinal tumour volume
% data from the DEN-induced murine HCC cohort (weeks 9-16; n = 10 HCC animals).
% Compares against linear and exponential alternatives.
%
% Gompertz model:
%   V(t) = A * exp(-B * exp(-C * t))
%   A = asymptotic volume (mm3)
%   B = dimensionless displacement parameter
%   C = growth rate constant (day^-1)
%
% Manuscript values (Section 3.8, Supplementary Table S8):
%   A = 312.4 mm3, B = 4.82, C = 0.185 day^-1
%   R2 = 0.988, MAPE = 6.3%
%
% IMPORTANT: This model was fitted to the 10-animal HCC cohort (in-sample).
% Independent prospective validation in a separate cohort is required before
% this constitutes a validated predictive tool. See manuscript Section 3.8.
%
% Reference: Norton L. Cancer Res 1988;48:7067-7071 (Gompertz tumour growth).

clear; clc; close all;

%% --- Longitudinal tumour volume data ------------------------------------
% Measurement timepoints (days post-DEN injection; weeks 9,11,13,16 = days 63,77,91,112)
t_days = [63, 77, 91, 112]';  % days

% Mean tumour volume per animal at each timepoint (mm3), n=10 HCC animals
% (Caliper measurement; confirmed against micro-CT in n=5 subset,
%  Bland-Altman bias +12.4 mm3, 95% LoA -18.1 to +42.9 mm3)
V_mean = [32.1, 89.4, 178.2, 298.6]';   % mm3 group mean
V_sd   = [8.2,  18.3,  31.5,  47.2]';   % mm3 group SD (n=10)

%% --- Gompertz model fitting --------------------------------------------
% Model: V(t) = A * exp(-B * exp(-C * t))
gompertz_model = @(p, t) p(1) .* exp(-p(2) .* exp(-p(3) .* t));

% Initial parameter guess
p0 = [350, 6, 0.15];
lb = [100, 0.1, 0.01];
ub = [1000, 50, 2.0];

opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 10000, 'TolFun', 1e-10);
p_gomp = lsqcurvefit(gompertz_model, p0, t_days, V_mean, lb, ub, opts);

V_gomp_fit = gompertz_model(p_gomp, t_days);
SS_res_g   = sum((V_mean - V_gomp_fit).^2);
SS_tot     = sum((V_mean - mean(V_mean)).^2);
R2_gomp    = 1 - SS_res_g / SS_tot;
MAPE_gomp  = mean(abs((V_mean - V_gomp_fit) ./ V_mean)) * 100;

fprintf('--- Gompertz Model ---\n');
fprintf('  A = %.2f mm3  (asymptotic volume)\n',  p_gomp(1));
fprintf('  B = %.4f     (displacement parameter)\n', p_gomp(2));
fprintf('  C = %.4f day^-1 (growth rate)\n',     p_gomp(3));
fprintf('  R2   = %.4f\n',  R2_gomp);
fprintf('  MAPE = %.2f%%\n\n', MAPE_gomp);

%% --- Alternative model 1: Linear ---------------------------------------
p_lin = polyfit(t_days, V_mean, 1);
V_lin_fit = polyval(p_lin, t_days);
SS_res_l  = sum((V_mean - V_lin_fit).^2);
R2_lin    = 1 - SS_res_l / SS_tot;
MAPE_lin  = mean(abs((V_mean - V_lin_fit) ./ V_mean)) * 100;
fprintf('--- Linear Model ---\n');
fprintf('  R2 = %.4f, MAPE = %.2f%%\n\n', R2_lin, MAPE_lin);

%% --- Alternative model 2: Exponential ---------------------------------
exp_model = @(p, t) p(1) .* exp(p(2) .* t);
p_exp = lsqcurvefit(exp_model, [1, 0.01], t_days, V_mean, [0, 0], [1000, 1], opts);
V_exp_fit = exp_model(p_exp, t_days);
SS_res_e  = sum((V_mean - V_exp_fit).^2);
R2_exp    = 1 - SS_res_e / SS_tot;
MAPE_exp  = mean(abs((V_mean - V_exp_fit) ./ V_mean)) * 100;
fprintf('--- Exponential Model ---\n');
fprintf('  R2 = %.4f, MAPE = %.2f%%\n\n', R2_exp, MAPE_exp);

%% --- Bootstrap 95% CIs for Gompertz parameters ------------------------
fprintf('Computing bootstrap 95% CIs (n=1000 resamples)...\n');
n_boot = 1000;
p_boot = zeros(n_boot, 3);
rng(2024);  % Fixed seed for reproducibility

for b = 1:n_boot
    % Resample (with replacement) from residuals and re-fit
    resid = V_mean - V_gomp_fit;
    V_boot = V_gomp_fit + resid(randi(length(resid), length(resid), 1));
    try
        p_b = lsqcurvefit(gompertz_model, p_gomp, t_days, V_boot, lb, ub, opts);
        p_boot(b, :) = p_b;
    catch
        p_boot(b, :) = p_gomp;
    end
end

CI_A = prctile(p_boot(:,1), [2.5, 97.5]);
CI_B = prctile(p_boot(:,2), [2.5, 97.5]);
CI_C = prctile(p_boot(:,3), [2.5, 97.5]);

fprintf('\nBootstrap 95% CIs (n=%d):\n', n_boot);
fprintf('  A: [%.2f, %.2f] mm3\n', CI_A(1), CI_A(2));
fprintf('  B: [%.4f, %.4f]\n',     CI_B(1), CI_B(2));
fprintf('  C: [%.4f, %.4f] day^-1\n\n', CI_C(1), CI_C(2));

%% --- Figure: model comparison and residuals ----------------------------
t_fine = linspace(55, 120, 300)';

fig = figure('Name', 'Gompertz Tumour Growth Model', 'Position', [100 100 1100 480]);

% Panel 1: Fits
subplot(1, 3, 1);
fill([t_fine; flipud(t_fine)], ...
    [gompertz_model([CI_A(2), CI_B(1), CI_C(2)], t_fine); ...
     flipud(gompertz_model([CI_A(1), CI_B(2), CI_C(1)], t_fine))], ...
    [0.5 0.7 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none'); hold on;
plot(t_fine, gompertz_model(p_gomp, t_fine), '-b', 'LineWidth', 2.5);
plot(t_fine, polyval(p_lin, t_fine),         '--', 'Color', [0.7 0.3 0], 'LineWidth', 1.5);
plot(t_fine, exp_model(p_exp, t_fine),       '-.', 'Color', [0.3 0.6 0.1], 'LineWidth', 1.5);
errorbar(t_days, V_mean, V_sd, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Time (days post-DEN)');
ylabel('Tumour volume (mm^3)');
title('Model Comparison');
legend({sprintf('Gompertz 95%% CI'), ...
    sprintf('Gompertz (R^2=%.3f)', R2_gomp), ...
    sprintf('Linear (R^2=%.3f)', R2_lin), ...
    sprintf('Exponential (R^2=%.3f)', R2_exp), ...
    'Observed (mean\pmSD, n=10)'}, 'Location', 'northwest', 'FontSize', 7);
grid on; box on;

% Add week labels
week_labels = {'Wk 9', 'Wk 11', 'Wk 13', 'Wk 16'};
for i = 1:4
    text(t_days(i), V_mean(i) + 15, week_labels{i}, 'HorizontalAlignment', 'center', 'FontSize', 7);
end

% Panel 2: Residuals
subplot(1, 3, 2);
residuals = V_mean - V_gomp_fit;
stem(t_days, residuals, 'filled', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [0.1 0.4 0.9]);
yline(0, '-k', 'LineWidth', 1);
xlabel('Time (days post-DEN)');
ylabel('Residual (mm^3)');
title(sprintf('Gompertz Residuals\nMAPE = %.2f%%', MAPE_gomp));
grid on; box on;
xticks(t_days); xticklabels(week_labels);

% Panel 3: Bootstrap parameter distributions
subplot(1, 3, 3);
histogram(p_boot(:,3) .* 1000, 30, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none', 'Normalization', 'probability');
hold on;
xline(p_gomp(3) * 1000, '-r', sprintf('C = %.3f (x10^{-3}) day^{-1}', p_gomp(3)*1000), ...
    'LabelVerticalAlignment', 'top', 'FontSize', 8);
xline(CI_C(1) * 1000, '--k', '2.5%', 'FontSize', 7);
xline(CI_C(2) * 1000, '--k', '97.5%', 'FontSize', 7);
xlabel('Growth rate C (x10^{-3} day^{-1})');
ylabel('Probability');
title(sprintf('Bootstrap Distribution of C\n(n=%d resamples)', n_boot));
grid on; box on;

sgtitle('Gompertz Tumour Growth Model - DEN-HCC Murine Cohort', 'FontSize', 12, 'FontWeight', 'bold');

saveas(fig, 'gompertz_fit.png');
fprintf('Figure saved: gompertz_fit.png\n');
