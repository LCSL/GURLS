function [param] = paramsel_fixsiglam_semi(X,y, opt)
% paramsel_fixlambda(X,y, OPT)
% Set the regularization parameter in OPT for semi-supervised learning
% 
% INPUTS:
% -OPT: not used
% 
% OUTPUT: struct with the following field:
% -lambdas: fixed value for the regularization parameter
% -sigma: fixed value for the kernel parameter
% -lambdas_m: fixed value for the manifold regularization parameter
% -sigma_m: fixed value for the manifold kernel parameter

%param.lambdas = 1e-3;
%param.sigma = 0.2;
%param.lambdas_m = 1;
%param.sigma_m = 0.1;

param.lambdas = 1e-5;
param.sigma = 0.3;
param.lambdas_m = 1e-3;
param.sigma_m = 0.1;
