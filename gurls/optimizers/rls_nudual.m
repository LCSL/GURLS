function [rls] = rls_nudual (X,y, opt)

% rls_nudual(X, y, opt)
% computes the regression function for the nu method in the dual space.
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

[n,T] = size(y);
nu=1;


if ~isfield(opt.paramsel,'Knorm');
    opt.paramsel.Knorm = norm(opt.kernel.K);
end

tau=1/(2*opt.paramsel.Knorm);

if isfield(opt.paramsel,'niter');
    niter = ceil(opt.paramsel.niter);
else
    niter = 0;
end

if isfield(opt.paramsel,'f0') && niter <= Niter;
    alpha1 = opt.paramsel.f0.alpha1;
    alpha2 = opt.paramsel.f0.alpha2;
else
    niter = 0;
    alpha2 = zeros(n,T);
    alpha1 = ((4*nu + 2)/(4*nu + 1)*tau)*y;
end


for i = niter:(Niter-1);
    u=((i-1)*(2*i-3)*(2*i+2*nu-1))/((i+2*nu-1)*(2*i+4*nu-1)*(2*i+2*nu-3));
    w=4*(((2*i+2*nu-1)*(i+nu-1)) /((i+2*nu-1)*(2*i+4*nu-1)) );
    alpha = alpha1 + u*(alpha1 - alpha2) +(w*tau)*(y - opt.kernel.K*alpha1);
    alpha2=alpha1;
    alpha1=alpha;
end

opt.paramsel.niter = Niter;
opt.paramsel.f0.alpha1 = alpha1;
opt.paramsel.f0.alpha2 = alpha2;

rls.C = alpha;
rls.X = X;
rls.W = [];
rls.CforInit = struct('Niter', Niter, 'alpha1',alpha1,'alpha2',alpha2);
