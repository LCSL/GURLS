function [cfr] = rls_smp(X,y,opt)

%	rls_smp(X,y,opt)
%	computes a classifier for the primal formulation of RLS.
%	The optimization is carried out using a stochastic mirror prox algorithm.
%	The regularization parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
%	NEEDS:
%		- opt.singlelamda
%		- opt.paramsel.lambdas
%		- opt.epochs
%		- opt.calibrate.gamma
%		- opt.Xte
%		- opt.yte


[n,d] = size(X);

T = size(y,2);

opt.cfr.W = zeros(d,T);
opt.cfr.R = zeros(d,T);
opt.cfr.W_sum = zeros(d,T);
opt.cfr.count = 0;
opt.cfr.acc_last = [];
opt.cfr.acc_avg = [];


% Run mulitple epochs
for i = 1:opt.epochs,
	opt.cfr = rls_smp_singlepass(X, y, opt);
end	

cfr = opt.cfr;
cfr.W = opt.cfr.W_sum/opt.cfr.count;
