function [cfr] = rls_primalrecinit(X, y, opt)
% rls_primalrecinit(X,y,opt)
% computes a classifier for the primal formulation of RLS. 
% The variables necessary for further recursive update are stored in the
% output structure.
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
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -Cinv: inverse of the regularized kernel matrix in the primal space
% -C: empty matrix
% -X: empty matrix

lambda = opt.singlelambda(opt.paramsel.lambdas);

if isfield(opt,'kernel');
    if isfield(opt.kernel,'XtX');
        XtX = opt.kernel.XtX;
    else
        XtX = X'*X;
    end
    if isfield(opt.kernel,'Xty');
        Xty = opt.kernel.Xty;
    else
        Xty = X'*y;
    end
    if isfield(opt.kernel,'n');
        n = opt.kernel.n;
    else
        n = size(y,1);
    end
else
    XtX = X'*X;
    Xty = X'*y;
    n = size(y,1);
end

d = size(XtX,1);

Cinv = pinv(XtX + (n*lambda)*eye(d));

cfr.W = Cinv*Xty;
cfr.C = [];
cfr.X = [];
cfr.Cinv = Cinv;
cfr.XtX = XtX;

