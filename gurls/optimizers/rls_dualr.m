function [cfr] = rls_dualr(X, y, opt)

%	rls_dualr(X,y,opt)
%	computes a classifier for the dual formulation of RLS.
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


n = size(opt.kernel.K,1);
T = size(y,2);

fprintf('\tSolving dual RLS...(n = %d, % = %d)', n, T);

K = opt.kernel.K + (n*lambda)*eye(n);

[Q,L,V] = tygert_svd(K,n);
Q = double(Q);
L = double(diag(L));

cfr.C = rls_eigen(Q, L, y, lambda,n);

if strcmp(opt.kernel.type, 'linear')
	cfr.W = X'*cfr.C;
	cfr.C = [];
	cfr.X = [];
else
	cfr.W = [];
end
