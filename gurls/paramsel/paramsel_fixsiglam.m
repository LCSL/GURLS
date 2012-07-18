function [param] = paramsel_fixsiglam(X, y, opt)
% paramsel_fixlambda(X,Y,OPT)
% Set the regularization parameter to the value set in OPT
% 
% INPUTS:
% -X: not used
% -Y: not used
% -OPT: not used
% 
% OUTPUT: struct with the following field:
% -lambdas: fixed value for the regularization parameter

param.lambdas = 10^-8;
param.sigma = 15;
