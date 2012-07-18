function [cfr] = rls_primalr(X, y, opt)
% rls_primalr(X,Y,OPT)
% computes a classifier for the primal formulation of RLS.
% It uses a randomized method to solve the associated linear system.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -X: input data matrix
% -Y: labels matrix
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: structure with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -C: empty matrix
% -X: empty matrix
% -dcomptime: time required for computing the singular value decomposition
% -rlseigtime: time required for computing the estimator given the SVD

lambda = opt.singlelambda(opt.paramsel.lambdas);

%fprintf('\tSolving primal RLS using Randomized SVD...\n');
[n,d] = size(X);

XtX = X'*X; % n\lambda is added in rls_eigen;

% tic;
[Q,L,U] = tygert_svd(XtX);
% cfr.dcomptime = toc;
Q = double(Q);
L = double(diag(L));

Xty = X'*y;

if isfield(opt,'W0') 
	Xty = Xty + opt.W0;
end

% tic;

cfr.W = rls_eigen(Q, L, Xty, lambda,d);
% cfr.rlseigtime = toc;

cfr.C = [];
cfr.X = [];
