%% === SETTINGS & MODEL ===
clear; close all; clc

% --- Input distributions (same as your original setup) ---
rw.param.dist.name   = "lognormal"; rw.param.mu.value    = 0.1;   rw.param.sigma.value = 0.02;
r.param.dist.name    = "uniform";   r.param.lower.value  = 50;    r.param.upper.value  = 150;
Tu.param.dist.name   = "uniform";   Tu.param.lower.value = 63070; Tu.param.upper.value = 115600;
Hu.param.dist.name   = "uniform";   Hu.param.lower.value = 990;   Hu.param.upper.value = 1110;
Tl.param.dist.name   = "uniform";   Tl.param.lower.value = 65;    Tl.param.upper.value = 1165;
Hl.param.dist.name   = "uniform";   Hl.param.lower.value = 700;   Hl.param.upper.value = 820;
L.param.dist.name    = "uniform";   L.param.lower.value  = 1120;  L.param.upper.value  = 1680;
Kw.param.dist.name   = "normal";    Kw.param.mu.value    = 10000; Kw.param.sigma.value = 1500;

% Sample size
nVar  = 8;
N     = 1e5;

% === (Important) Your sortcop function must be in the path ===
% function [Xs, Cs]=sortcop(C,X)
%   ... (your implementation) ...

%% === INDEPENDENT LHS (SUPS) AND MAPPING TO MARGINALS ===
SUPS = lhsdesign(N, nVar);     % stratified in [0,1]

SRVS = zeros(N, nVar);         % real-value samples (independent)
SRVS(:,1) = icdf(rw.param.dist.name, SUPS(:,1), rw.param.mu.value, rw.param.sigma.value);
SRVS(:,2) = icdf(r.param.dist.name,   SUPS(:,2), r.param.lower.value, r.param.upper.value);
SRVS(:,3) = icdf(Tu.param.dist.name,  SUPS(:,3), Tu.param.lower.value, Tu.param.upper.value);
SRVS(:,4) = icdf(Hu.param.dist.name,  SUPS(:,4), Hu.param.lower.value, Hu.param.upper.value); % Hu
SRVS(:,5) = icdf(Tl.param.dist.name,  SUPS(:,5), Tl.param.lower.value, Tl.param.upper.value);
SRVS(:,6) = icdf(Hl.param.dist.name,  SUPS(:,6), Hl.param.lower.value, Hl.param.upper.value); % Hl
SRVS(:,7) = icdf(L.param.dist.name,   SUPS(:,7), L.param.lower.value,  L.param.upper.value);
SRVS(:,8) = icdf(Kw.param.dist.name,  SUPS(:,8), Kw.param.mu.value,    Kw.param.sigma.value);

%% === CORRELATION MATRIX FOR THE COPULA (Hu–Hl only) ===
% Variable order: {rw, r, Tu, Hu, Tl, Hl, L, Kw}
C = eye(nVar);
rho = 0.7;
C(4,6) = rho; C(6,4) = rho;   % correlation only between Hu (4) and Hl (6)

%% === IMPOSE DEPENDENCE VIA GAUSSIAN COPULA (sortcop) ===
[Xcorr, Cs] = sortcop(C, SRVS);     % Xcorr keeps marginals, imposes rank correlation
LHS.SRVS = Xcorr;                   % use the correlated set from here on

% Quick check (optional):
disp('Sample correlation (Pearson) after copula, Cs ='); disp(Cs([4 6],[4 6]));

%% === MODEL EVALUATION ===
% Assuming you have borehole.m vectorized:
y = borehole(LHS.SRVS(:,1), LHS.SRVS(:,2), LHS.SRVS(:,3), LHS.SRVS(:,4), ...
             LHS.SRVS(:,5), LHS.SRVS(:,6), LHS.SRVS(:,7), LHS.SRVS(:,8));

%% === SUMMARY STATISTICS ===
mu_hat       = mean(y);
sigma_hat    = std(y);
var_hat      = var(y);
mc_se        = sigma_hat/sqrt(N);

fprintf('Output f: mean=%.6g, std=%.6g, var=%.6g, MC-SE=%.6g\n', mu_hat, sigma_hat, var_hat, mc_se);

%% === DIAGNOSTICS OF MARGINALS AND SCATTERS (correlated) ===
varNames = {'rw','r','Tu','Hu','Tl','Hl','L','Kw'};
X = LHS.SRVS;

% Marginal histograms
figure('Name','Marginal Histograms — Correlated Inputs','Color','w');
tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
for j=1:nVar
    nexttile; histogram(X(:,j),'Normalization','pdf');
    title(varNames{j}); xlabel(varNames{j}); ylabel('PDF (emp.)'); grid on
end

% Scatter matrix (reduced sample for visualization)
idx = randperm(N, min(3000,N));
Xplot = X(idx,:);
figure('Name','Scatter Matrix — Correlated Inputs','Color','w');
plotmatrix(Xplot);
sgtitle('LHS + Copula (Hu–Hl with \rho=0.7)');

%% === SENSITIVITY ANALYSIS: Random Sampling + Binning (Scatter + binned curve) ===
% Idea: for each variable X_j, plot scatter of (X_j, y) and overlay
% the mean of y per bin in X_j (approximate E[y|X_j]).

nBins = 20;             % number of bins for the line
maxPtsScatter = 8000;   % limit scatter points for readability

figure('Name','Sensitivity — Random Sampling + Binning','Color','w');
tiledlayout(2,4,'Padding','compact','TileSpacing','compact');

for j = 1:nVar
    nexttile;
    % Subsample for scatter:
    ii = randperm(N, min(maxPtsScatter, N));
    scatter(X(ii,j), y(ii), 6, '.', 'MarkerEdgeAlpha',0.25); hold on; grid on
    
    % Bins: use quantile edges (robust to scale/outliers)
    qEdges = quantile(X(:,j), linspace(0,1,nBins+1));
    % Avoid duplicate edges (e.g. for discrete vars)
    qEdges = unique(qEdges);
    if numel(qEdges) < 3
        % fallback: linspace bins
        qEdges = linspace(min(X(:,j)), max(X(:,j)), nBins+1);
    end
    % Bin centers
    binCenters = 0.5*(qEdges(1:end-1)+qEdges(2:end));
    yBinMean   = nan(size(binCenters));
    for b = 1:numel(binCenters)
        inb = X(:,j) >= qEdges(b) & X(:,j) < qEdges(b+1);
        if b == numel(binCenters)
            % include last right edge
            inb = X(:,j) >= qEdges(b) & X(:,j) <= qEdges(b+1);
        end
        if any(inb)
            yBinMean(b) = mean(y(inb));
        end
    end
    plot(binCenters, yBinMean, 'LineWidth', 2);   % binned trend E[y|X_j]
    xlabel(varNames{j}); ylabel('f_{borehole}');
    title([varNames{j} ' vs f(x) — binning']);
    legend('Random sampling','Binned mean','Location','best'); legend boxoff
    hold off
end
sgtitle('Sensitivity via Random Sampling + Binning (with Hu–Hl dependence)');

