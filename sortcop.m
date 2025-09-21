function [Xs, Cs]=sortcop(C,X)
% generating rank correlated samples using Gaussian copula
% Input:
%   C    : correlation matrix of the variables (nvar,nvar)
%   X : uncorrelated samples
% Output:
%   Xs     : Rank correlated samples with C
%   Cs: Correlatoin matrix of rank-sorted samples of Xs.
%   Gurkan Sin (2016), DTU

[nsample nvar]=size(X);
u = copularnd('Gaussian',C,nsample);
for i=1:nvar
    [s,ii(:,i)] = sort(u(:,i));
    Xs(ii(:,i),i) = sort(X(:,i));
end
Cs=corr(Xs);
end
