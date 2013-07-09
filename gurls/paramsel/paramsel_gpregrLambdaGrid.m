function vout = paramsel_gpregrLambdaGrid(X,y,opt)
% paramsel_gpregrLambdaGrid(X,Y,OPT)
% Performs parameter selection for gaussian process regression by 
% maximizing the likelihood.
%
% INPUTS:
% -X: input data matrix
% -Y: labels matrix
% -OPT: structure of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- kernel.K (set by the kernel_* routines)
%   fields with default values set through the defopt function:
%       - nlambda
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -guesses: array of guesses for the regularization parameter lambda 
% -perf: matrix with the likelihood for each value of the regularization
%        parameter and for each class
% -lambdas: array of values of the regularization parameter lambda
%           minimizing the validation error for each class

if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
end

n = size(y,1);
        
if isfield(opt,'lambdamin')
    lmin = opt.lambdamin;
else
    lmin = 0.001;
end
if isfield(opt,'lambdamax')
    lmax = opt.lambdamax;
else
    lmax = 10;
end
powers = linspace(0,1,tot);
guesses = lmin.*(lmax/lmin).^(powers);
	

% L = sort(L,'descend');
% % maximum eigenvalue
% lmax = L(1);
% CumSumEig = cumsum(L);
% firstPercentile = find(CumSumEig./CumSumEig(end)>.999,1,'first');
% lmin = max(L(firstPercentile), 200*sqrt(eps));
% powers = linspace(0,1,tot);
% guesses = lmin.*(lmax/lmin).^(powers);

for i = 1:tot
    opt.paramsel.lambdas = guesses(i);
    opt.rls = rls_gpregr(X,y,opt);
    perf(i,:) = opt.rls.logprob;
    vout.datafit(i,:) = opt.rls.datafit;
    vout.penalty(i,:) = opt.rls.penalty;

end
[dummy,idx] = max(perf,[],1);

vout.lambdas = guesses(idx);
vout.perf = perf;
vout.guesses = guesses;

