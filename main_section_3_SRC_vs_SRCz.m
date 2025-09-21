%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Borehole Function — Uncertainty & SRC Sensitivity (single-file version)
%
% Workflow:
%   1) Define input distributions
%   2) LHS in unit probability space (SUPS) + inverse CDF to real space (SRVS)
%   3) Evaluate model (borehole)
%   4) Linear regression on standardized variables -> SRCs + R^2
%   5) Sorted table + bar plot (saved for LaTeX)
%
% NOTE (lognormal): icdf('Lognormal',p,mu,sigma) usa os parâmetros da
%   Normal subjacente de log(X). Se rw.mu/sigma estiverem no domínio REAL
%   (média e desvio de rw), converta para o domínio LOG (veja trecho abaixo).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
rng(1,'twister');  % reprodutibilidade

%% 1) Input uncertainty definition

% Borehole radius (lognormal: mu,sigma in LOG-domain)
rw.param.dist.name   = "lognormal";
rw.param.mu.value    = 0.1;   % ASSUMED log-domain mu
rw.param.sigma.value = 0.02;  % ASSUMED log-domain sigma

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

% Variable order (must match borehole)
varNames = {'rw','r','Tu','Hu','Tl','Hl','L','Kw'};
m = numel(varNames);

% Sample size
N = 1e5;

%% 2) LHS sampling: SUPS -> SRVS (inverse CDF mapping)

% (A) Unit probability space (N x m)
SUPS = lhsdesign(N, m);

% (B) Real value space
SRVS = zeros(N, m);

% --- If rw.mu/sigma are REAL-domain (mean m_rw, std s_rw), convert to LOG-domain:
% m_rw = rw.param.mu.value; s_rw = rw.param.sigma.value;
% sigma_log = sqrt(log(1 + (s_rw/m_rw)^2));
% mu_log    = log(m_rw) - 0.5*sigma_log^2;
% SRVS(:,1) = icdf('lognormal', SUPS(:,1), mu_log, sigma_log);

% Using provided (assumed log-domain) parameters:
SRVS(:,1) = icdf(rw.param.dist.name, SUPS(:,1), rw.param.mu.value, rw.param.sigma.value);
SRVS(:,2) = icdf(r.param.dist.name,  SUPS(:,2), r.param.lower.value,  r.param.upper.value);
SRVS(:,3) = icdf(Tu.param.dist.name, SUPS(:,3), Tu.param.lower.value, Tu.param.upper.value);
SRVS(:,4) = icdf(Hu.param.dist.name, SUPS(:,4), Hu.param.lower.value, Hu.param.upper.value);
SRVS(:,5) = icdf(Tl.param.dist.name, SUPS(:,5), Tl.param.lower.value, Tl.param.upper.value);
SRVS(:,6) = icdf(Hl.param.dist.name, SUPS(:,6), Hl.param.lower.value, Hl.param.upper.value);
SRVS(:,7) = icdf(L.param.dist.name,  SUPS(:,7), L.param.lower.value,  L.param.upper.value);
SRVS(:,8) = icdf(Kw.param.dist.name, SUPS(:,8), Kw.param.mu.value,    Kw.param.sigma.value);

% Alias columns for readability (rw, r, Tu, Hu, Tl, Hl, L, Kw)
X = SRVS;
rw_s = X(:,1); r_s  = X(:,2); Tu_s = X(:,3); Hu_s = X(:,4);
Tl_s = X(:,5); Hl_s = X(:,6); L_s  = X(:,7); Kw_s = X(:,8);

%% 3) Monte Carlo model evaluations (vectorized)

y = borehole(rw_s, r_s, Tu_s, Hu_s, Tl_s, Hl_s, L_s, Kw_s);

% Summary statistics
mu_hat    = mean(y);
sigma_hat = std(y);
var_hat   = var(y);
mc_err    = sigma_hat / sqrt(N);
fprintf('Mean=%.6g, Std=%.6g, Var=%.6g, MCerr=%.6g (N=%d)\n', mu_hat, sigma_hat, var_hat, mc_err, N);

%% 4) SRCs via standardized regression (no intercept) — RECOMMENDED

% Standardize X and y (z-score)
[Xz,~,~] = zscore(X);
[yz,~,~] = zscore(y);

% On standardized variables, no intercept -> coefficients are SRCs
SRC  = Xz \ yz;                 % m x 1
yhat = Xz * SRC;
R2z  = corr(yhat, yz)^2;

% Order by |SRC|
[~,idxSort]   = sort(abs(SRC), 'descend');
inputs_sorted = varNames(idxSort);
SRC_sorted    = SRC(idxSort);

% Table & print
T = table(varNames(:), SRC, 'VariableNames', {'Input','SRC'});
disp('--- SRC results (standardized, no intercept) ---');
disp(T);
fprintf('R^2 (standardized, no intercept) = %.4f\n', R2z);

%% 5) Plot SRCs (ordered) — ProSafe/DTU color scheme

dtu_blue  = [0, 58,108]/255;
dark_gray = [74,74,74]/255;

figure('Name','SRCs (Standardized Regression Coefficients)','Color','w');
cat_inputs = categorical(inputs_sorted);
cat_inputs = reordercats(cat_inputs, inputs_sorted);
bar(cat_inputs, SRC_sorted, 'FaceColor', dtu_blue);
ylabel('SRC'); xlabel('Input');
title(sprintf('Standardized Regression Coefficients (R^2 = %.3f)', R2z));
grid on;
exportgraphics(gcf, 'src_barplot.png','Resolution',300); % For LaTeX

%% (Optional) Alternative check — unstandardized + intercept (do NOT mix in report)

B     = [ones(N,1) X] \ y;           % intercept + slopes
b     = B(2:end);
yhat0 = [ones(N,1) X]*B;
R2    = corr(yhat0, y)^2;
SRC_u = b .* (std(X,0,1)'/std(y));   % classic SRC definition

fprintf('R^2 (unstandardized + intercept) = %.4f; ||SRC_u - SRC||_2 = %.3g\n', ...
        R2, norm(SRC_u - SRC, 2));


[br bint] = regress(y,[ones(N,1) X]);

ymgs = [ones(N,1) X] * br ;
Rrgs = corr(ymgs,y);
Rr2gs = Rrgs^2 
br(1)=[]; % delete intercept
SRCr=br.*std(X)'./std(y)


sum(SRCr.^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 4: Linear regression - A linear model is fitted to MC simulation results (ym) and input data (X) to obtain the Standardized Regression Coefficients (SRCs) for sensitivity analysis
% Settings for lsqlin function 
lb=ones(m,1)* -1 ;% SRC = [-1,1]
ub=ones(m,1)*  1 ;%
options = optimoptions(@lsqlin, 'MaxIterations', 3e7);

% Obtaining coefficients b using linear least squares regression with bounds
[b,resnorm,residual,exitflag,output] = lsqlin(X,y,[],[],[],[],lb,ub,[],options);


% Checking LSQ convergence
disp('LSQ Convergence')
disp(exitflag)

% New linear model predictions with regression coefficients
ym = X * b; 

% Checking goodness of fit 
R = corr(ym,y);
R2 = R^2 ;
disp('Goodness of fit (R_squared)')
disp(R2) % R2 < 0.7 => the model canNOT be linearized sufficiently => SRCs are NOT reliable sensitivity measures => use Spearman's rank correlation or Sobol indices instead

% Evaluating SRC and sum of squared SRCs
SRC = b.*std(X)'./std(y)
unity = sum(SRC.^2) % Sum of SCR_squared should be equal to 1 for linear models

