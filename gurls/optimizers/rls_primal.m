function [cfr] = rls_primal (X, y, opt)

% rls_primal(X,y,opt)
% computes a classifier for the primal formulation of RLS.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
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


indices = 1:size(X,1);
if isprop(opt,'split_fixed_indices') && isprop(opt,'notTrainOnValidation') && opt.notTrainOnValidation
    indices = opt.split_fixed_indices;
end

n = numel(indices);

% check if matrices XtX and Xty have been previously computed during
% parameter selection
if isfield(opt.paramsel,'XtX') && n == size(X,1);
    XtX = opt.paramsel.XtX;
else
    XtX = X(indices,:)'*X(indices,:); % d x d matrix.
end
if isfield(opt.paramsel,'XtX') && n == size(X,1);
    Xty = opt.paramsel.Xty;
else
    Xty = X(indices,:)'*y(indices,:); % d x T matrix.
end
cfr.W = rls_primal_driver( XtX, Xty, n, lambda );
cfr.C = [];
cfr.X = [];

