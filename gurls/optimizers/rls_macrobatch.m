function [cfr] = rls_macrobatch(X, Y, opt)
% PRIMAL
% Assumes d-reasonably small, but n-very large
% Implements a warm restart iterative macrobatch 
% descent optimization.
% Each independent classifier is solved in closed form,
% using rls_primal or dual.


opt.W0 = 0;
nIter = 20;
batchsize = 6000; % multiple of 'd'?

[n,d] = size(X);

idx = unique(randsample(n,batchsize));

opt.expname = 'bagging_dummy';
opt.savefile = ['bagging_dummy.mat'];

try
	opt.seq = opt.bagseq;
	opt.process = opt.bagproc;
catch
	error('Specify opt.bagseq and opt.bagproc for training individual models\n');
end

opt.singlelambda = @median;

for i = 1:nIter,
	fprintf('\n\t\t[macrobatch %d/%d]', i,nIter);

	% Subsample - x, 
	sub_X = X(idx,:);
	sub_Y = Y(idx,:);

	gurls(sub_X, sub_Y, opt,1);

	load(opt.savefile);
	opt.W0 = opt.rls.W;
end
