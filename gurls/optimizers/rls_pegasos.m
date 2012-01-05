function [cfr] = rls_pegasos(X, bY, opt)

%	rls_pegasos(X,y,opt)
%	computes a classifier for the primal formulation of RLS.
%	The optimization is carried out using a stochastic gradient descent algorithm.
%	The regularization parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
%	NEEDS:
%		- opt.singlelamda
%		- opt.paramsel.lambdas
%		- opt.epochs
%		- opt.Xte
%		- opt.yte
	
% Pegasos-type algorithm for RLS
% Solves RLS in the primal
% X : n x d matrix data
% bY = binary coded Y
% opt = options

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
