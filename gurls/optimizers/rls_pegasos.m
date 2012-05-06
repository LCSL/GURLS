function [cfr] = rls_pegasos(X, bY, opt)
% rls_pegasos(X,BY,OPT)
% computes a classifier for the primal formulation of RLS.
% The optimization is carried out using a stochastic gradient descent algorithm.
% The regularization parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -X: input data matrix
% -BY: binary coded labels matrix
% -OPT: structure of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%       - epochs
%   fields that need to be added by hand
%       -Xte
%       -yte
%   fields with default values set through the defopt function:
%		- singlelambda
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: structure with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -t0: stepsize parameter 
% -W_sum: sum of the classifiers across iterations
% -count: number of iterations
% -acc_last: accuracy of the solution computed in the last iteration
% -acc_avg: average accuracy across iterations

[n,d] = size(X);

T = size(bY,2);

opt.cfr.W = zeros(d,T);
opt.cfr.W_sum = zeros(d,T);
opt.cfr.count = 0;
opt.cfr.acc_last = [];
opt.cfr.acc_avg = [];


% Run mulitple epochs
for i = 1:opt.epochs,
	if opt.cfr.count == 0
		opt.cfr.t0 = ceil(norm(X(1,:))/sqrt(opt.singlelambda(opt.paramsel.lambdas)));
		fprintf('\n\tt0 is set to : %f\n', opt.cfr.t0);
	end
	opt.cfr = rls_pegasos_singlepass(X, bY, opt);
end	
cfr = opt.cfr;
cfr.W = opt.cfr.W_sum/opt.cfr.count;
