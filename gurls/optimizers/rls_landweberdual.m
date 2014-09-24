function [rls] = rls_landweberdual(X,y, opt)

% rls_landweberdual(X, y, opt)
% computes the regression function for landweber regularization in the dual space.
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

if isfield(opt.paramsel,'niter');
    niter = ceil(opt.paramsel.niter);
else
    niter = 0;
end

if isfield(opt.paramsel,'f0') && niter <= Niter;
    alpha = opt.paramsel.f0;
else
    niter = 0;
    alpha=zeros(n,T);
end

if ~isfield(opt.paramsel,'Knorm');
    opt.paramsel.Knorm = norm(opt.kernel.K); 
end

tau=1/(2*opt.paramsel.Knorm); 


for i = niter:(Niter-1);
    alpha = alpha + tau*(y- opt.kernel.K*alpha);
end

opt.paramsel.f0 = alpha;
opt.paramsel.niter = Niter;

rls.C = alpha;
rls.X = X;
rls.W = [];



