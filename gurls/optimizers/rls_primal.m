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

n = size(y,1);

% check if matrices XtX and Xty have been previously computed during
% parameter selection
if isfield(opt.paramsel,'XtX');
    XtX = opt.paramsel.XtX;
else
    XtX = X'*X; % d x d matrix.
end
if isfield(opt.paramsel,'XtX');
    Xty = opt.paramsel.Xty;
else
    Xty = X'*y; % d x T matrix.
end
cfr.W = rls_primal_driver( XtX, Xty, n, lambda );
cfr.C = [];
cfr.X = [];

