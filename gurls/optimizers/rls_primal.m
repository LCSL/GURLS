function [cfr] = rls_primal (X, y, opt)

%	rls_primal(X,y,opt)
%	computes a classifier for the primal formulation of RLS.
%	The regularization parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
%	NEEDS:
%		- opt.singlelamda
%		- opt.paramsel.lambdas

lambda = opt.singlelambda(opt.paramsel.lambdas);

fprintf('\tSolving primal RLS...\n');

[n,d] = size(X);

XtX = X'*X; % d x d matrix.
Xty = X'*y; % d x T matrix.

cfr.W = rls_primal_driver( XtX, Xty, n, lambda );
cfr.C = [];
cfr.X = [];

