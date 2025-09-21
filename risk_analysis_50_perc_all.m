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

%% 3) Specify input distributions (reduced intervals by 50%)
% Borehole function input parameters and their distributions
% Each variable has fields: distribution type + parameters (mu, sigma, bounds)

% Borehole radius (lognormal) → halve sigma
rw.param.dist.name   = "lognormal";
rw.param.mu.value    = 0.1;
rw.param.sigma.value = 0.5 * 0.02;   % 50% reduction

% Radius of influence (uniform) → halve interval around midpoint
r.param.dist.name    = "uniform";
mid = (50+150)/2; range = (150-50)/2; newRange = 0.5*range;
r.param.lower.value  = mid - newRange;
r.param.upper.value  = mid + newRange;

% Transmissivity of upper aquifer (uniform)
Tu.param.dist.name   = "uniform";
mid = (63070+115600)/2; range = (115600-63070)/2; newRange = 0.5*range;
Tu.param.lower.value = mid - newRange;
Tu.param.upper.value = mid + newRange;

% Potentiometric head of upper aquifer (uniform)
Hu.param.dist.name   = "uniform";
mid = (990+1110)/2; range = (1110-990)/2; newRange = 0.5*range;
Hu.param.lower.value = mid - newRange;
Hu.param.upper.value = mid + newRange;

% Transmissivity of lower aquifer (uniform)
Tl.param.dist.name   = "uniform";
mid = (65+1165)/2; range = (1165-65)/2; newRange = 0.5*range;
Tl.param.lower.value = mid - newRange;
Tl.param.upper.value = mid + newRange;

% Potentiometric head of lower aquifer (uniform)
Hl.param.dist.name   = "uniform";
mid = (700+820)/2; range = (820-700)/2; newRange = 0.5*range;
Hl.param.lower.value = mid - newRange;
Hl.param.upper.value = mid + newRange;

% Length of borehole (uniform)
L.param.dist.name    = "uniform";
mid = (1120+1680)/2; range = (1680-1120)/2; newRange = 0.5*range;
L.param.lower.value  = mid - newRange;
L.param.upper.value  = mid + newRange;

% Hydraulic conductivity of borehole (normal) → halve sigma
Kw.param.dist.name   = "normal";
Kw.param.mu.value    = 10000;
Kw.param.sigma.value = 0.5 * 1500;   % 50% reduction



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
% Given: y is the N-by-1 vector of model outputs (m^3/yr) from Task 1
thr = 9500;                   % critical flow rate (m^3/yr)
N   = numel(fborehole.param.values);

p_exceed = mean(fborehole.param.values > thr);     % Monte Carlo estimate
% 95% binomial CI for exceedance probability
z = 1.96;
se = sqrt(max(p_exceed*(1-p_exceed)/N, 0));
ci_baseline = [p_exceed - z*se, p_exceed + z*se];

fprintf('Baseline P(f>%.0f) = %.4f (95%% CI: [%.4f, %.4f])\n', ...
        thr, p_exceed, ci_baseline(1), ci_baseline(2));

 %%   7.1 Risk Quantification 

 %fborehole.param.values
thr = 9500;
pd  = fitdist(fborehole.param.values , 'Normal');        % estima mu e sigma
p_exceed_norm = 1 - cdf(pd, thr);  % igual a 1 - normcdf(thr, pd.mu, pd.sigma)
p_exceed_norm 


