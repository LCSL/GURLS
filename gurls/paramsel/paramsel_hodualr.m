function [vout] = paramsel_hodualr(X, y, opt)

%	paramsel_hodualr(X,y,opt)
%	Performs parameter selection when the dual formulation of RLS is used.
%	The eigendecomposition used to compute the regularization path is computed using a randoamized method.
%	The performance measure specified by opt.hoperf is maximized.
%
%	NEEDS:	
%		- opt.split
%		- opt.nholdouts
%		- opt.kernel.type
%		- opt.kernel.K
%		- opt.nlambda
%		- opt.hoperf

savevars = [];
for nh = 1:opt.nholdouts
	if strcmp(class(opt.split),'cell')
		tr = opt.split{nh}.tr;
		va = opt.split{nh}.va;
	else	
		tr = opt.split.tr;
		va = opt.split.va;
	end	


	[n,T]  = size(y(tr,:));
	
	if strcmp(opt.kernel.type,'load')
		d = n; % Here U is same as Q
	else
		[n, d] = size(X(tr,:));
	end
	
	[Q,L,U] = tygert_svd(opt.kernel.K(tr,tr),n);
	Q = double(Q);
	L = double(diag(L));
	
	guesses = paramsel_lambdaguesses(L, min(n,d), n, opt);
	tot = opt.nlambda;
	ap = zeros(tot,T);
	
	QtY = Q'*y(tr,:);
	for i = 1:tot
		%guesses(i) = lmin*(q^i);
		%%%%%% REPLICATING CODE FROM RLS_DUAL %%%%%%%%
		opt.rls.C = rls_eigen(Q,L,QtY,guesses(i),n);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		if size(X,1) > 0
			Xva = X(va,:);
			opt.rls.X = X(tr,:);
		else
			Xva = [];
		end
	
		yva = y(va,:);
		if strcmp(opt.kernel.type,'linear')
			opt.rls.W = X(tr,:)'*opt.rls.C; 
		elseif strcmp(opt.kernel.type,'load')
			opt.predkernel.type = 'load';
			opt.predkernel.K = opt.kernel.K(va,tr);
		else
			opt.predkernel = predkernel_traintest(X(va,:),y(va,:),opt);
		end	
	
		opt.pred = pred_dual(Xva,yva,opt);
		opt.perf = opt.hoperf(Xva,yva,opt);
		%p{i} = perf(scores,yho,{'precrec'});
		for t = 1:T
			ap(i,t) = opt.perf.forho(t);
		end	
	end	
	[dummy,idx] = max(ap,[],1);	
	vout.lambdas_round{nh} = guesses(idx);
	vout.forho{nh} = ap;
	vout.guesses{nh} = guesses;
	% This is awesome
	if numel(savevars) > 0
		[ST,I] = dbstack();
		save(ST(1).name,savevars{:});
	end	
end	

if numel(vout.lambdas_round) > 1
	lambdas = cell2mat(vout.lambdas_round');
	vout.lambdas = mean(lambdas);
else
	vout.lambdas = vout.lambdas_round{1};
end
