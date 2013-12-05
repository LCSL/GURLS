function [rls] = rls_landweberprimal (X, y, opt)

% rls_landweberprimal(X,y,opt)
% computes the regression function for landweber regularization in the primal space.
% The regularization parameter (i.e. the number of iterations) is set to the one found in opt.paramsel.
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
% -C: empty matrix
% -X: empty matrix



Niter = opt.singlelambda(opt.paramsel.lambdas);

d = size(X,2);
T = size(y,2);

if isfield(opt.paramsel,'f0');
    W = opt.paramsel.f0;
else
    W=zeros(d,T);
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

for i = 1:Niter;
    W = W + tau*(Xty- XtX*W);
end

rls.C = [];
rls.X = [];
rls.W = W;

