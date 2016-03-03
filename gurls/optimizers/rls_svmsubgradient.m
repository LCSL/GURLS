function [rls] = rls_svmsubgradient(X,y, opt)

% rls_svmsubgradient(X, y, opt)
% computes the regression function for svm subgradient method.
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

lambda = opt.singlelambda(opt.paramsel.lambdas);
Niter = opt.Niter;

[n,T] = size(y);

if isfield(opt.paramsel,'niter');
    niter = max(1,ceil(opt.paramsel.niter));
else
    niter = 1;
end

if isfield(opt.paramsel,'f0') && niter <= Niter;
    alpha = opt.paramsel.f0;
else
    niter = 1;
    alpha=zeros(n,T);
end

iter = niter:(Niter);
gamma = 1./iter;

for i = iter;
    yhat = opt.kernel.K*alpha;
    yxyhat = y.*yhat;
    indic = (yxyhat <= 1);
    alpha = alpha*(1-2*gamma(i)*lambda) + (gamma(i)/n)*(y.*indic);
end

opt.paramsel.f0 = alpha;
opt.paramsel.niter = Niter;

if strcmp(opt.kernel.type, 'linear')
	rls.W = X'*alpha;
	rls.C = [];
	rls.X = [];
else
	rls.W = [];
    rls.C = alpha;
	rls.X = X;
end



