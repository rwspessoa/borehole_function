%% Workflow
% WORKFLOW — UNCERTAINTY ANALYSIS (Monte Carlo + LHS) FOR THE BOREHOLE FUNCTION
%
% 1) Define the model:
%    - State the borehole function f(x) and check units/consistency.
%
% 2) Implement the model:
%    - Provide a vectorized function (e.g., borehole.m) callable with arrays.
%
% 3) Specify input distributions:
%    - Inputs: rw, r, Tu, Hu, Tl, Hl, L, Kw.
%    - Assign probability models (Uniform/Normal/Lognormal) with parameters/bounds.
%
% 4) Choose sample size(s):
%    - Select N for Latin Hypercube Sampling (e.g., 1e3, 1e4, 1e5) to study convergence.
%
% 5) Generate LHS samples:
%    - Produce stratified quantiles per input and map via inverse CDF to each distribution.
%
% 6) Run model evaluations:
%    - Evaluate f for all N samples; store outputs in a vector y (size N x 1).
%
% 7) Compute summary stats:
%    - Mean (mu_hat), Std (sigma_hat), Variance (sigma2_hat = sigma_hat^2).
%
% 8) Plot distribution diagnostics:
%    - Histogram of y; empirical CDF of y.
%
% 9) Convergence analysis:
%    - Running mean/variance vs N;
%    - Monte Carlo error ~ 1/sqrt(N);
%    - Mean/Std (M/S) vs N.
%
% 10) Interpretation:
%    - Comment on adequacy of N, stability of estimates, and salient features of y.


clear  
close all  
clc   

%% 1) Define the model 

%% 2) Implement the model 

    % code: borehole.m

%% 3) Specify input distributions
% Borehole function input parameters and their distributions
% Each variable has fields: distribution type + parameters (mu, sigma, bounds)

% Borehole radius
rw.param.dist.name   = "lognormal";
rw.param.mu.value    = 0.1;
rw.param.sigma.value = 0.02;

% Radius of influence
r.param.dist.name    = "uniform";
r.param.lower.value  = 50;
r.param.upper.value  = 150;

% Transmissivity of upper aquifer
Tu.param.dist.name   = "uniform";
Tu.param.lower.value = 63070;
Tu.param.upper.value = 115600;

% Potentiometric head of upper aquifer
Hu.param.dist.name   = "uniform";
Hu.param.lower.value = 990;
Hu.param.upper.value = 1110;

% Transmissivity of lower aquifer
Tl.param.dist.name   = "uniform";
Tl.param.lower.value = 65;
Tl.param.upper.value = 1165;

% Potentiometric head of lower aquifer
Hl.param.dist.name   = "uniform";
Hl.param.lower.value = 700;
Hl.param.upper.value = 820;

% Length of borehole
L.param.dist.name    = "uniform";
L.param.lower.value  = 1120;
L.param.upper.value  = 1680;

% Hydraulic conductivity of borehole
Kw.param.dist.name   = "normal";
Kw.param.mu.value    = 10000;
Kw.param.sigma.value = 1500;



%% 4) Choose sample size(s):

LHS.param.variable.size = 8;
LHS.param.variable.sample = 1e5;


%% 5) Generate LHS samples:
      % Sampling in Unit Probability Space
        LHS.SUPS = lhsdesign(LHS.param.variable.sample, ...
                             LHS.param.variable.size);

      % Sampling in Real Value Space
        LHS.SRVS = zeros(LHS.param.variable.sample, ...
                         LHS.param.variable.size);




      % project probability space into real-valued input 
      % space using distribution-specific inverse transforms

      LHS.SRVS(:,1) = icdf(rw.param.dist.name, ...
                           LHS.SUPS(:,1), ...
                           rw.param.mu.value, ...
                           rw.param.sigma.value);

      LHS.SRVS(:,2) = icdf(r.param.dist.name, ...
                           LHS.SUPS(:,2), ...
                           r.param.lower.value, ...
                           r.param.upper.value);

      LHS.SRVS(:,3) = icdf(Tu.param.dist.name, ...
                           LHS.SUPS(:,3), ...
                           Tu.param.lower.value, ...
                           Tu.param.upper.value);

      LHS.SRVS(:,4) = icdf(Hu.param.dist.name, ...
                           LHS.SUPS(:,4), ...
                           Hu.param.lower.value, ...
                           Hu.param.upper.value);

      LHS.SRVS(:,5) = icdf(Tl.param.dist.name, ...
                           LHS.SUPS(:,5), ...
                           Tl.param.lower.value, ...
                           Tl.param.upper.value);

      LHS.SRVS(:,6) = icdf(Hl.param.dist.name, ...
                           LHS.SUPS(:,6), ...
                           Hl.param.lower.value, ...
                           Hl.param.upper.value);

      LHS.SRVS(:,7) = icdf(L.param.dist.name, ...
                           LHS.SUPS(:,7), ...
                           L.param.lower.value, ...
                           L.param.upper.value);

      LHS.SRVS(:,8) = icdf(Kw.param.dist.name, ...
                           LHS.SUPS(:,8), ...
                           Kw.param.mu.value, ...
                           Kw.param.sigma.value);



% ==== Grid de histogramas e CDF ====
X = LHS.SRVS;
varNames = {'rw','r','Tu','Hu','Tl','Hl','L','Kw'};
nVar = numel(varNames);

% Histogramas
figure('Name','Marginal Histograms — LHS.SRVS','Color','w');
tiledlayout(2,4, 'Padding','compact', 'TileSpacing','compact');
for j = 1:nVar
    nexttile;
    histogram(X(:,j), 'Normalization','pdf');
    title(varNames{j});
    xlabel(varNames{j}); ylabel('PDF (emp.)'); grid on;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define ProSafe/DTU color scheme
dtu_red     = [144, 26, 30] / 255;
dtu_blue    = [0, 58, 108] / 255;
light_gray  = [242, 242, 242] / 255;
dark_gray   = [74, 74, 74] / 255;

% Sample for plotting
N = size(LHS.SRVS,1);
idx = randperm(N, min(3000, N));
Xplot = LHS.SRVS(idx,:);

% Variables names
varNames = {'rw','r','Tu','Hu','Tl','Hl','L','Kw'};

% Create scatter plot matrix
figure('Name','Scatter Plot Matrix - LHS.SRVS','Color','w');
[h, ax] = plotmatrix(Xplot);

% Customize scatter plot colors and marker style
for i = 1:numel(h)
    set(h(i), 'Marker', '.', 'MarkerEdgeColor', dtu_blue, 'MarkerSize', 8);
end

% Customize diagonal histograms
% These are the ax on diagonal – histogram axes are not directly returned in h
% We use ax handles: ax(i,i) is the histogram plot for var i
nVar = length(varNames);
for i = 1:nVar
    ax_ii = ax(i, i);
    hold(ax_ii, 'on');
    hHist = histogram(ax_ii, Xplot(:, i), 'FaceColor', dtu_blue, 'EdgeColor', dark_gray);
    % overlay original histogram – you might want to clear previous bins
    % or adjust nearly blank diagonal.
    ax_ii.Color = light_gray;  % background color
    hold(ax_ii, 'off');
end

sgtitle('Scatter Plot Matrix of LHS Samples','FontWeight','Bold','FontSize',14);

% Label axes
for k = 1:length(varNames)
    xlabel(ax(end, k), varNames{k}, 'Color', dark_gray, 'FontSize',12);
    ylabel(ax(k, 1), varNames{k}, 'Color', dark_gray, 'FontSize',12);
end


%% 6) Run model evaluations:

fborehole.param.values = borehole(LHS.SRVS(:,1), ...
                                  LHS.SRVS(:,2), ... 
                                  LHS.SRVS(:,3), ...
                                  LHS.SRVS(:,4), ...
                                  LHS.SRVS(:,5), ...
                                  LHS.SRVS(:,6), ... 
                                  LHS.SRVS(:,7), ...
                                  LHS.SRVS(:,8));

%% 7) Compute summary stats:

fborehole.param.mu = mean(fborehole.param.values);
fborehole.param.sigma = std(fborehole.param.values);
fborehole.param.var = var(fborehole.param.values);
fborehole.param.mcstderror = std(fborehole.param.values)/sqrt(length(fborehole.param.values));



%% 8) Plot distribution diagnostics:

% ==== Scatter plots: each input variable vs. model output ====
varNames = {'rw','r','Tu','Hu','Tl','Hl','L','Kw'};
X = LHS.SRVS;
y = fborehole.param.values;

figure('Name','Inputs vs Output','Color','w');
tiledlayout(2,4, 'Padding','compact', 'TileSpacing','compact');

for j = 1:numel(varNames)
    nexttile;
    scatter(X(1:5000,j), y(1:5000), 8, [144 26 30]/255, '.'); % DTU red
    xlabel(varNames{j});
    ylabel('f_{borehole}');
    title([varNames{j} ' vs f(x)']);
    grid on;
end
sgtitle('Scatter Plots — Inputs vs Borehole Function Output');




% ==== Histogram of borehole function output ====
y = fborehole.param.values;

figure('Name','Histogram of Borehole Output','Color','w');
histogram(y, 50, ...                          % 50 bins
    'FaceColor', [144 26 30]/255, ...         % DTU red
    'EdgeColor', 'w', ...
    'Normalization','pdf');                   % normalizado para PDF

xlabel('f_{borehole} (m^3/yr)');
ylabel('Probability Density');
title('Distribution of Borehole Function Output');
grid on;

hold on;
pd = fitdist(y,'Normal');
x_vals = linspace(min(y), max(y), 200);
plot(x_vals, pdf(pd,x_vals), 'k-', 'LineWidth', 1.5);
legend('Empirical Histogram','Normal Fit');
hold off;

exportgraphics(gcf, 'histogram_fborehole.png', 'Resolution',300);






%% 9) Convergence analysis (with saving figures for LaTeX)

% A)  without restart LHS 
% B)  with restart LHS

%--------------------------------------------------------------------------

% % % A)  without restart LHS 
y = fborehole.param.values(:);    % ensure column
N  = numel(y);
n  = (1:N)';                      % running sample sizes

% cumulative sums
S1 = cumsum(y);
S2 = cumsum(y.^2);

mu_run = S1 ./ n;                               
var_pop_est = S2./n - mu_run.^2;                    
var_unb = var_pop_est .* (n ./ max(n-1,1));         
sigma_run = sqrt(max(var_unb, 0));                 
mc_se = sigma_run ./ sqrt(n);

z = 1.96;
ci_lo = mu_run - z .* mc_se;
ci_hi = mu_run + z .* mc_se;

MS_ratio = mu_run ./ sigma_run;

k = unique(round(logspace(1, log10(N), 400)));     
kk = k(k >= 2);                                    

% --- Figure 1: Running mean + 95% CI ---
f1 = figure('Name','Convergence: Running Mean with 95% CI','Color','w');
plot(n(kk), mu_run(kk), 'LineWidth', 1.4); hold on;
plot(n(kk), ci_lo(kk), '--', n(kk), ci_hi(kk), '--', 'LineWidth', 1.0);
xlabel('Sample size N'); ylabel('Running mean \mu(N)');
title('Running Mean of f_{borehole} with 95% Confidence Band');
grid on; legend('Mean','95% CI bounds','Location','best');
exportgraphics(f1, 'convergence_running_mean_ci.png','Resolution',300);

% --- Figure 2: Running variance ---
f2 = figure('Name','Convergence: Running Variance','Color','w');
plot(n(kk), var_unb(kk), 'LineWidth', 1.4);
xlabel('Sample size N'); ylabel('Running variance \sigma^2(N)');
title('Running Variance of f_{borehole}');
grid on;
exportgraphics(f2, 'convergence_running_variance.png','Resolution',300);

% --- Figure 3: MC error ~ N^{-1/2} ---
f3 = figure('Name','Convergence: Monte Carlo SE vs N','Color','w');
loglog(n(kk), mc_se(kk), 'LineWidth', 1.4); hold on;
c_ref = mc_se(kk(end)) * sqrt(n(kk(end)));
ref_line = c_ref ./ sqrt(n(kk));
loglog(n(kk), ref_line, '--', 'LineWidth', 1.2);
xlabel('N (log scale)'); ylabel('MC standard error \sigma/\sqrt{N}');
title('Monte Carlo Standard Error and N^{-1/2} Reference');
grid on; legend('MC SE','\propto N^{-1/2}','Location','best');
exportgraphics(f3, 'convergence_mc_se_loglog.png','Resolution',300);

% --- Figure 4: Mean/Std (M/S) vs N ---
f4 = figure('Name','Convergence: M/S vs N','Color','w');
plot(n(kk), MS_ratio(kk), 'LineWidth', 1.4);
xlabel('Sample size N'); ylabel('M/S = \mu(N) / \sigma(N)');
title('Mean-to-Std Ratio vs Sample Size');
grid on;
exportgraphics(f4, 'convergence_ms_ratio.png','Resolution',300);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % B)  with restart LHS 
% ============================
% 2) Define sample sizes
% ============================
Nmax = 1e5;
sample_sizes = 2.^(5:ceil(log2(Nmax))); % 2^10=1024 up to >=1e5

% Preallocate results
mu_vals    = zeros(size(sample_sizes));
sigma_vals = zeros(size(sample_sizes));
var_vals   = zeros(size(sample_sizes));
mc_err     = zeros(size(sample_sizes));

% ============================
% 3) Loop over sample sizes
% ============================
for k = 1:numel(sample_sizes)
    N = sample_sizes(k);

    % Latin Hypercube Sampling in probability space
    SUPS = lhsdesign(N, 8);

    % Map to real value space
    SRVS = zeros(N,8);
    SRVS(:,1) = icdf(rw.param.dist.name, SUPS(:,1), rw.param.mu.value, rw.param.sigma.value);
    SRVS(:,2) = icdf(r.param.dist.name, SUPS(:,2), r.param.lower.value, r.param.upper.value);
    SRVS(:,3) = icdf(Tu.param.dist.name, SUPS(:,3), Tu.param.lower.value, Tu.param.upper.value);
    SRVS(:,4) = icdf(Hu.param.dist.name, SUPS(:,4), Hu.param.lower.value, Hu.param.upper.value);
    SRVS(:,5) = icdf(Tl.param.dist.name, SUPS(:,5), Tl.param.lower.value, Tl.param.upper.value);
    SRVS(:,6) = icdf(Hl.param.dist.name, SUPS(:,6), Hl.param.lower.value, Hl.param.upper.value);
    SRVS(:,7) = icdf(L.param.dist.name, SUPS(:,7), L.param.lower.value, L.param.upper.value);
    SRVS(:,8) = icdf(Kw.param.dist.name, SUPS(:,8), Kw.param.mu.value, Kw.param.sigma.value);

    % Evaluate model
    y = borehole(SRVS(:,1), SRVS(:,2), SRVS(:,3), SRVS(:,4), ...
                 SRVS(:,5), SRVS(:,6), SRVS(:,7), SRVS(:,8));

    % Statistics
    mu_vals(k)    = mean(y);
    sigma_vals(k) = std(y);
    var_vals(k)   = var(y);
    mc_err(k)     = sigma_vals(k) / sqrt(N);

    fprintf('N=%d: mean=%.4f, std=%.4f, MCerr=%.4f\n', ...
        N, mu_vals(k), sigma_vals(k), mc_err(k));
end

% ============================
% 4) Plot convergence curves
% ============================
f5 = figure('Color','w');
subplot(2,2,1);
semilogx(sample_sizes, mu_vals,'-o','LineWidth',1.2);
xlabel('Sample size N (log)'); ylabel('Mean'); grid on;
title('Convergence of mean');

subplot(2,2,2);
semilogx(sample_sizes, sigma_vals,'-o','LineWidth',1.2);
xlabel('Sample size N (log)'); ylabel('Std dev'); grid on;
title('Convergence of std dev');

subplot(2,2,3);
loglog(sample_sizes, mc_err,'-o','LineWidth',1.2,'DisplayName','MC error'); hold on;

% Reference line ~ N^{-1/2}
ref = mc_err(end) * sqrt(sample_sizes(end)) ./ sqrt(sample_sizes);
loglog(sample_sizes, ref, '--','LineWidth',1.5,'Color',[0 58 108]/255,'DisplayName','N^{-1/2}');

xlabel('Sample size N (log scale)');
ylabel('MC error');
title('Monte Carlo standard error');
grid on;
legend('show','Location','northeast');

subplot(2,2,4);
semilogx(sample_sizes, mu_vals ./ sigma_vals,'-o','LineWidth',1.2);
xlabel('Sample size N (log)'); ylabel('M/S ratio'); grid on;
title('Mean-to-Std ratio');
exportgraphics(f5,'convergence_restart_LHS.png','Resolution',300);
