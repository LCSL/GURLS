function [rls] = rls_landweberprimal_lgs (X,y,opt)

% rls_landweberprimal(X, y, opt)
% computes the regression function for landweber regularization in the primal space.
% The regularization parameter (i.e. the number of iterations) is set to the one found in opt.paramsel.
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

Niter = ceil(opt.singlelambda(opt.paramsel.lambdas));

d = size(X,2);
T = size(y,2);

if isfield(opt.paramsel,'niter');
    niter = ceil(opt.paramsel.niter);
else
    niter = 1;
end

if isfield(opt.paramsel,'f0') && niter <= Niter;
    W = opt.paramsel.f0;
else
    W = zeros(d,T);
end

if isfield(opt.paramsel,'XtX');
    XtX = opt.paramsel.XtX;
else
    XtX = X'*X; % d x d matrix.
end
if isfield(opt.paramsel,'Xty');
    Xty = opt.paramsel.Xty;
else
    Xty = X'*y; % d x T matrix.
end

if isfield(opt.paramsel,'XtXnorm');
    tau=1/(2*opt.paramsel.XtXnorm); 
else
    tau=1/(2*norm(XtX));
end

for i = niter:Niter;
    W = W + tau*(X'*(y./(1+exp(y.*(X*W)))));
end

opt.paramsel.f0 = W;
opt.paramsel.niter = Niter;

rls.C = [];
rls.X = [];
rls.W = W;