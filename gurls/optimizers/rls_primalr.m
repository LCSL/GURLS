function [cfr] = rls_primalr(X, y, opt)

%	rls_primalr(X,y,opt)
%	computes a classifier for the primal formulation of RLS.
% 	It uses a randomized method to solve the associated linear system.
%	The regularization parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
%	NEEDS:
%		- opt.singlelamda
%		- opt.paramsel.lambdas
%		- opt.kernel.K
%		- opt.kernel.type

lambda = opt.singlelambda(opt.paramsel.lambdas);

fprintf('\tSolving primal RLS using Randomized SVD...\n');
[n,d] = size(X);

XtX = X'*X; % n\lambda is added in rls_eigen;

tic;
[Q,L,U] = tygert_svd(XtX,d);
cfr.dcomptime = toc;
Q = double(Q);
L = double(diag(L));

Xty = X'*y;

if isfield(opt,'W0') 
	Xty = Xty + opt.W0;
end

tic;

cfr.W = rls_eigen(Q, L, Xty, lambda,d);
cfr.rlseigtime = toc;

cfr.C = [];
cfr.X = [];
