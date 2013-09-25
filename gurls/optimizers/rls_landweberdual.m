function [rls] = rls_landweberdual(X, y, opt)

% rls_landweberdual(X,y,opt)
% computes the regression function for landweber regularization in the dual space.
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

[n,T] = size(y);

if isfield(opt.paramsel,'f0');
    alpha = opt.paramsel.f0;
else
    alpha=zeros(n,T);
end

if isfield(opt.paramsel,'Knorm');
    tau=1/(2*opt.paramsel.Knorm); 
else
    tau=1/(2*norm(opt.kernel.K)); 
end

for i = Niter;
    alpha = alpha + tau*(y- opt.kernel.K*alpha);
end

rls.C = alpha;
rls.X = X;
rls.W = [];



