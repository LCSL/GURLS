function [out] = rls_lapfilter(X,y,opt)
		if strcmp(class(opt.split),'cell')
			tr = opt.split{nh}.tr;
			va = opt.split{nh}.va;
		else	
			tr = opt.split.tr;
			va = opt.split.va;
		end	

		[n,T] = size(y(opt.split.tr,:));
		opt.rls.C = zeros(n,T);
		opt.rls.X = X(opt.split.tr,:);

		oldperf = zeros(1,T);
		topScoring = zeros(n,T);
		epsilon = opt.singlelambda(opt.paramsel.lambdas);
		epsilon = 1;

		J = eye(numel(opt.split.tr));
		opt.kernel.K = opt.kernel.K + epsilon*J;
		opt.laplacian.L = opt.laplacian.L + epsilon*J;

		opt.predkernel = predkernel_traintest(X(va,:),y(va,:),opt);
		
		y(opt.split.u,:) = 0;
		ytr = y(opt.split.tr,:);
		J(numel(opt.split.l)+1:end,:) = 0;

		LK = opt.laplacian.L*opt.kernel.K;
		LK = pinv(LK);

		for i = 1:opt.maxIter
			opt.rls.C = opt.rls.C + LK*(ytr - J*opt.kernel.K*opt.rls.C);
		%	opt.rls.C = opt.rls.C + J*K(ytr-J*opt.kernel.K*opt.rls.C) + 
			opt.pred = pred_dual(X(va,:),y(va,:),opt);
			opt.perf = opt.hoperf(X(va,:),y(va,:),opt);
			for t = 1:T
				if opt.perf.forho(t) > oldperf(t)
					topScoring(:,t) = opt.rls.C(:,t);
					bestIter(t) = i;
					oldperf = opt.perf.forho;
				end
			end	
		end
		out = opt.rls;
		out.bestIter = bestIter;
end		
		
				
