function [cfr] = rls_dual (X, y, opt)
% rls_dual(X, y, opt) computes a classifier for the dual formulation of RLS.
%
% The regularization parameter is set to the one found in opt.paramsel. In
% case of multiclass problems, the regularizers need to be combined with
% the opt.singlelambda function.
%
% INPUTS:
% -  X: input data matrix
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
% OUTPUT: struct with the following fields:
% if kernel.type='linear'
% -W: matrix of coefficient vectors of primal rls estimator for each class
% -C: empty matrix
% -X: empty matrix
% else
% -W: empty matrix
% -C: matrix of coefficient vectors of dual rls estimator for each class
% -X: input data matrix

lambda = opt.singlelambda(opt.paramsel.lambdas);

indices = 1:size(X,1);
if isprop(opt,'split_fixed_indices') && isprop(opt,'notTrainOnValidation') && opt.notTrainOnValidation
    indices = opt.split_fixed_indices;
end

n = numel(indices);

%fprintf('\tSolving dual RLS...(n = %d, % = %d)', n, T);

try
    R = chol(opt.kernel.K(indices,indices) + (n*lambda)*eye(n));
    
    cfr.C = R\(R'\y(indices,1));
catch
    [Q,L,~] = svd(opt.kernel.K(indices,indices));
    Q = double(Q);
    L = double(diag(L));
    cfr.C = rls_eigen(Q, L, Q'*y(indices,:), lambda, n);
end

if strcmp(opt.kernel.type, 'linear')
    cfr.W = X(indices,:)'*cfr.C;
    cfr.C = [];
    cfr.X = [];
else
    cfr.W = [];
    cfr.X = X(indices,:);
end
