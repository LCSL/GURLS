function [out] = rls_monline(X,y,opt)

[n,T] = size(y);
d = size(X,2);
opt.rls.W0 = zeros(d,T);

oldperf = zeros(1,T);
topScoring = zeros(d,T);

bestIter = ones(1,T);

for i = 1:opt.nSplits

	opt.paramsel.lambdas = (T/sqrt(i)) * opt.paramsel.lambdas;

	opt.rls = rls_monline_driver(X(opt.split.monline{i}.idx,:), y(opt.split.monline{i}.idx,:), opt);

	opt.pred = pred_primal(X(opt.split.va,:),y(opt.split.va,:),opt);
	opt.perf = opt.hoperf(X(opt.split.va,:),y(opt.split.va,:),opt);

	for t = 1:T
		if opt.perf.forho(t) > oldperf(t)
			topScoring(:,t) = opt.rls.W(:,t);
			bestIter(t) = i;
			oldperf = opt.perf.forho;
		end
	end	
	
	opt.rls.W0 = opt.rls.W;

end	

out.W = topScoring;
out.bestIter = bestIter;
