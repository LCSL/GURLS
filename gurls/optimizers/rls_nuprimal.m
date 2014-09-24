function [rls] = rls_nuprimal (X,y, opt)
% rls_nuprimal(X, y, opt)
% computes the regression function for the nu method in the primal space.
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
nu=1;

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

if isfield(opt.paramsel,'niter');
    niter = ceil(opt.paramsel.niter);
else
    niter = 1;
end

if isfield(opt.paramsel,'f0') && niter <= Niter;
    beta1 = opt.paramsel.f0.beta1;
    beta2 = opt.paramsel.f0.beta2;
else
    beta2=zeros(d,T);
    beta1 = ((4*nu + 2)/(4*nu + 1)*tau)*Xty;
end


for i = 1:Niter;
    beta = beta1 + u*(beta1 - beta2) +(w*tau)*(Xty - XtX*beta1);
    u=((i-1)*(2*i-3)*(2*i+2*nu-1))/((i+2*nu-1)*(2*i+4*nu-1)*(2*i+2*nu-3));
    w=4*(((2*i+2*nu-1)*(i+nu-1)) /((i+2*nu-1)*(2*i+4*nu-1)) );
    beta2=beta1;
    beta1=beta;
end

opt.paramsel.f0.beta1 = beta1;
opt.paramsel.f0.beta2 = beta2;

opt.paramsel.niter = Niter;

rls.W = beta;
rls.X = [];
rls.C = [];
rls.WforInit = struct('beta1',beta1,'beta2',beta2);
