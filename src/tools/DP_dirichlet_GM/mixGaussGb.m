function [label, Theta, w, llh] = mixGaussGb(X, opt)
% Collapsed Gibbs sampling for Dirichlet process (infinite) Gaussian mixture model (a.k.a. DPGM).
% This is a wrapper function which calls underlying Dirichlet process mixture model.
% Input:
%   X: d x n data matrix
%   opt(optional): prior parameters
% Output:
%   label: 1 x n cluster label
%   Theta: 1 x k structure of trained Gaussian components
%   w: 1 x k component weight vector
%   llh: loglikelihood
% Written by Mo Chen (sth4nth@gmail.com).

if isfield(opt,'seed')
    rng(opt.seed);
else
    rng default
end
[d,n] = size(X);
mu = mean(X,2);
Xo = bsxfun(@minus,X,mu);
s = sum(Xo(:).^2)/(d*n);
if nargin == 1
    kappa0 = 1;
    m0 = mean(X,2);
    nu0 = d;
    S0 = s*eye(d);
    alpha0 = 1;
else
    if isfield(opt,'kappa')
        kappa0 = opt.kappa;
    else
        kappa0 = 1;
    end
    if isfield(opt,'m')
        m0 = opt.m;
    else
        m0 = mean(X,2);
    end
    if isfield(opt,'nu')
        nu0 = opt.nu;
    else
        nu0 = d;
    end
    if isfield(opt,'S')
        S0 = opt.S;
    else
        S0 = s*eye(d);
    end
    if isfield(opt,'alpha')
        alpha0 = opt.alpha;
    else
        alpha0 = 1;
    end
end
prior = GaussWishart(kappa0,m0,nu0,S0);
[label, Theta, w, llh] = mixDpGb(X,alpha0,prior);
