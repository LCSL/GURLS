function [cfr] = rls_dualr(X,y, opt)
% rls_dualr(X, y, opt)
% computes a classifier for the dual formulation of RLS.
% It uses a randomized method to solve the associated linear system.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -OPT: struct of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routine)
%       - kernel (set by the kernel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
%		- kernel.type
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: structure with the following fields:
% if kernel.type='linear'
% -W: matrix of coefficient vectors of primal rls estimator for each class
% -C: empty matrix
% -X: empty matrix
% else
% -W: empty matrix
% -C: matrix of coefficient vectors of dual rls estimator for each class
% -X: empty matrix

lambda = opt.singlelambda(opt.paramsel.lambdas);


n = size(opt.kernel.K,1);
T = size(y,2);

%fprintf('\tSolving dual RLS...(n = %d, % = %d)', n, T);

K = opt.kernel.K;

k = max(1,round(opt.eig_percentage*n/100));
[Q,L,V] = tygert_svd(K,k);
Q = double(Q);
L = double(diag(L));

cfr.C = rls_eigen(Q, L, Q'*y, lambda,n);

if strcmp(opt.kernel.type, 'linear')
	cfr.W = X'*cfr.C;
	cfr.C = [];
	cfr.X = [];
else
	cfr.W = [];
	cfr.X = X;
end
