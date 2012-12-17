function [cfr] = rls_randfeats (X, y, opt)

% rls_primalRF(X,y,opt)
% computes a classifier for the primal formulation of RLS using the Random
% Features approach proposed in: 
%   Ali Rahimi, Ben Recht;
%   Random Features for Large-Scale Kernel Machines;
%   in Neural Information Processing Systems (NIPS) 2007.
%
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
%       - randfeat.D
%       - randfeat.samplesize
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -C: empty matrix
% -X: empty matrix

lambda = opt.singlelambda(opt.paramsel.lambdas);

%fprintf('\tSolving primal RLS...\n');

[n,d] = size(X);

cfr.W = rls_primal_driver( opt.kernel.XtX, opt.kernel.Xty, n, lambda );

cfr.C = [];
cfr.X = [];

