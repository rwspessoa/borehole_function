%% Sobol' sensitivity analysis (Saltelli scheme) for the Borehole function
clear; close all; clc;
rng(1,'twister');   % reproducibility for non-quasi steps

% ========== 1) Define input distributions ==========
% Order must match borehole.m: {rw, r, Tu, Hu, Tl, Hl, L, Kw}

rw.param.dist.name   = "lognormal"; 
rw.param.mu.value    = 0.1;   % assumed LOG-domain mu
rw.param.sigma.value = 0.02;  % assumed LOG-domain sigma

r.param.dist.name    = "uniform";   
r.param.lower.value  = 50;    
r.param.upper.value  = 150;

Tu.param.dist.name   = "uniform";   
Tu.param.lower.value = 63070; 
Tu.param.upper.value = 115600;

Hu.param.dist.name   = "uniform";   
Hu.param.lower.value = 990;   
Hu.param.upper.value = 1110;

Tl.param.dist.name   = "uniform";   
Tl.param.lower.value = 65;    
Tl.param.upper.value = 1165;

Hl.param.dist.name   = "uniform";   
Hl.param.lower.value = 700;   
Hl.param.upper.value = 820;

L.param.dist.name    = "uniform";   
L.param.lower.value  = 1120;  
L.param.upper.value  = 1680;

Kw.param.dist.name   = "normal";    
Kw.param.mu.value    = 10000; 
Kw.param.sigma.value = 1500;

varNames = {'rw','r','Tu','Hu','Tl','Hl','L','Kw'};
m = numel(varNames);

% ========== 2) Sampling parameters ==========
% N base samples; total evaluations ≈ (m+2)N
% Choose N as a power of 2 for Sobol sequences
N = 1e5; 
ss = sobolset(m,'Skip',1e3,'Leap',1e2); 
ss = scramble(ss,'MatousekAffineOwen');

% Base matrices A and B in [0,1]
A_u = net(ss, N);                     
ss2 = sobolset(m,'Skip',5e4,'Leap',1e2); 
ss2 = scramble(ss2,'MatousekAffineOwen');
B_u = net(ss2, N);

% ========== 3) Utility function: map [0,1] -> real space ==========
map_to_real = @(U) [ ...
    icdf(rw.param.dist.name, U(:,1), rw.param.mu.value,    rw.param.sigma.value), ...
    icdf(r.param.dist.name,  U(:,2), r.param.lower.value,  r.param.upper.value), ...
    icdf(Tu.param.dist.name, U(:,3), Tu.param.lower.value, Tu.param.upper.value), ...
    icdf(Hu.param.dist.name, U(:,4), Hu.param.lower.value, Hu.param.upper.value), ...
    icdf(Tl.param.dist.name, U(:,5), Tl.param.lower.value, Tl.param.upper.value), ...
    icdf(Hl.param.dist.name, U(:,6), Hl.param.lower.value, Hl.param.upper.value), ...
    icdf(L.param.dist.name,  U(:,7), L.param.lower.value,  L.param.upper.value), ...
    icdf(Kw.param.dist.name, U(:,8), Kw.param.mu.value,    Kw.param.sigma.value) ];

% NOTE: if rw.mu/sigma are REAL-domain values, convert to LOG-domain first
% m_rw = 0.1; s_rw = 0.02;
% sigma_log = sqrt(log(1 + (s_rw/m_rw)^2));
% mu_log    = log(m_rw) - 0.5*sigma_log^2;
% SRVS(:,1) = icdf('lognormal', U(:,1), mu_log, sigma_log);

% ========== 4) Build mixed matrices A_B^(i) ==========
A = A_u;  
B = B_u;            

Y_A   = zeros(N,1);
Y_B   = zeros(N,1);
Y_ABi = zeros(N,m);

% Evaluate A and B (after mapping to real space)
XA = map_to_real(A);
XB = map_to_real(B);

Y_A = borehole(XA(:,1), XA(:,2), XA(:,3), XA(:,4), XA(:,5), XA(:,6), XA(:,7), XA(:,8));
Y_B = borehole(XB(:,1), XB(:,2), XB(:,3), XB(:,4), XB(:,5), XB(:,6), XB(:,7), XB(:,8));

% For each i, build A_B^(i): all columns from A except column i, which is from B
for i = 1:m
    ABi = A;                     
    ABi(:,i) = B(:,i);           
    XABi = map_to_real(ABi);     
    Y_ABi(:,i) = borehole(XABi(:,1), XABi(:,2), XABi(:,3), XABi(:,4), ...
                          XABi(:,5), XABi(:,6), XABi(:,7), XABi(:,8));
end

% ========== 5) Sobol' estimators ==========
VarY = var([Y_A; Y_B], 1);   % population variance for stability

% First-order indices (Saltelli 2010)
S_first = zeros(m,1);
for i = 1:m
    S_first(i) = mean( Y_B .* (Y_ABi(:,i) - Y_A) ) / VarY;
end

% Total-order indices (Jansen 1999 / Saltelli 2010)
S_total = zeros(m,1);
for i = 1:m
    S_total(i) = mean( (Y_A - Y_ABi(:,i)).^2 ) / (2*VarY);
end

% ========== 6) Reporting and plotting ==========
% Sort by importance (total effects)
[~, idxT] = sort(S_total, 'descend');
varNames_sortedT = varNames(idxT);
S_first_sorted   = S_first(idxT);
S_total_sorted   = S_total(idxT);

% ProSafe/DTU colors
dtu_blue  = [0, 58,108]/255;
dtu_red   = [144,26, 30]/255;

% Bar chart of first-order vs total indices
figure('Name','Sobol Indices','Color','w');
bar(categorical(varNames_sortedT), [S_first_sorted, S_total_sorted], 'grouped');
legend('First-order S_i','Total S_{Ti}','Location','best');
ylabel('Sensitivity index');
xlabel('Input');
title(sprintf('Sobol\\'' Indices (N=%d base, total evals ≈ %d)', N, (m+2)*N));
grid on;
exportgraphics(gcf, 'sobol_indices.png','Resolution',300);

% Print results as table
T_sobol = table(varNames', S_first, S_total, 'VariableNames', {'Input','S_first','S_total'});
disp('--- Sobol'' indices (unsorted) ---'); 
disp(T_sobol);

% ========== 7) (Optional) Compare with SRC to check linearity ==========
X = XA;                          
y = Y_A;
[Xz,~,~] = zscore(X); [yz,~,~] = zscore(y);
SRC  = Xz \ yz;                  
R2z  = corr(Xz*SRC, yz)^2;

fprintf('\nLinearity check: R^2_z (SRC model) = %.3f\n', R2z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local function: borehole
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function f = borehole(rw, r, Tu, Hu, Tl, Hl, L, Kw)
%    % Borehole function: flow rate through a borehole
%    log_term = log(r ./ rw);
%    bracket  = 1 + (2 .* L .* Tu) ./ (log_term .* rw.^2 .* Kw) + (Tu ./ Tl);
%    f = (2 .* pi .* Tu .* (Hu - Hl)) ./ (log_term .* bracket);
%end
