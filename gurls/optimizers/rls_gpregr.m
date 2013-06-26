function [cfr] = rls_gpregr(X, y, opt)
% rls_gpregr(X,y,opt)
% Performs gaussian process regression
% 
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: struct of options with the following fields (and subfields):
%       -opt.paramsel.lambdas (set through the paramsel routines)
%       -opt.singlelambda (default)
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
%  -L
%  -alpha
%  -X
%  -logprob

noise = opt.singlelambda(opt.paramsel.lambdas);

n = size(opt.kernel.K,1);

cfr.L = chol(opt.kernel.K + noise^2*eye(n));

cfr.alpha = cfr.L\(cfr.L'\y);
cfr.X = X;
for t = 1:size(y,2);
    cfr.logprob(t) = -(y(:,t)'*cfr.alpha(:,t) + 2*sum(log(diag(cfr.L))) + n*log(2*pi))/2;
    cfr.datafit(t) = -(y(:,t)'*cfr.alpha(:,t))/2;
    cfr.penalty(t) = -sum(log(diag(cfr.L)));
end

