function [vout] = paramsel_hoprimal(X,y,opt)

%	paramsel_hoprimal(X,y,opt)
%	Performs parameter selection when the primal formulation of RLS is used.
%	The performance measure specified by opt.hoperf is maximized.
%
%	NEEDS:	
%		- opt.split
%		- opt.nholdouts
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

	if opt.hoMOnline
		tr = opt.split.monline{1}.idx;
	end	



	[n,d] = size(X(tr,:));
	[n,T]  = size(y(tr,:));
	
	K = X(tr,:)'*X(tr,:);
	
	tot = opt.nlambda;
	[Q,L] = eig(K);
	Q = double(Q);
	L = double(diag(L));
	QtXtY = Q'*(X(tr,:)'*y(tr,:));
	
	guesses = paramsel_lambdaguesses(L, min(n,d), n, opt);
	
	% Replaced with new function paramsel_lambdaguesses
	%filtered = L(L > 200*eps^0.5);
	%lmin = min(filtered)/n;
	%lmax = max(filtered)/n;
	%q = (lmax/lmin)^(1/tot);
	%guesses = zeros(1,tot);
	
	ap = zeros(tot,T);
	for i = 1:tot
		% guesses(i) = lmin*(q^i); % replaced the computation with new function paramsel_lambdaguesses
		opt.rls.W = rls_eigen(Q,L,QtXtY,guesses(i),n);
		opt.pred = pred_primal(X(va,:),y(va,:),opt);
		opt.perf = opt.hoperf(X(va,:),y(va,:),opt);
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
