function [vout] = paramsel_hoprimal(X,y,opt)

%	bigparamsel_dhoprimal(X,y,opt)
%	Performs parameter selection when the primal formulation of RLS is used.
%	This specific implementation assumes XtXm, Xty, XvayXva and Xvatyva
%	have already been computed using gdm and are available in files.
%	The performance measure specified by opt.hoperf is maximized.
%
%	NEEDS:	
%		- opt.nlambda
%		- opt.hoperf
%		- opt.files.Xva_filename
%		- opt.files.yva_filename
%		- opt.files.XtX_filename
%		- opt.files.Xty_filename
%		- opt.files.XvatXva_filename
%		- opt.files.Xvatyva_filename

	X.Transpose(true);
	y.Transpose(true);
	
	%% Load pre-computed products
	
	Xva = bigarray.Obj(opt.files.Xva_filename);
	yva = bigarray.Obj(opt.files.yva_filename);
	
	Xva.Transpose(true);
	yva.Transpose(true);
	
	t = load(opt.files.XtX_filename);	XtX = t.data;
	t = load(opt.files.Xty_filename);	Xty = t.data;
	t = load(opt.files.XvatXva_filename); 	XvatXva = t.data;
	t = load(opt.files.Xvatyva_filename); 	Xvatyva = t.data;;
	clear t;
	
	%%
	
	K = XtX - XvatXva;
	Xty = Xty - Xvatyva;
	
	n = X.NumItems - Xva.NumItems;
	d = Sizes(X,1);
	T = Sizes(y,1);
	
	
	tot = opt.nlambda;
	[Q,L] = eig(K);
	Q = double(Q);
	L = double(diag(L));
	QtXtY = Q'*Xty;
	
	guesses = paramsel_lambdaguesses(L, min(n,d), n, opt);
	
	
	ap = zeros(tot,T);
	for i = 1:tot
		opt.rls.W = rls_eigen(Q,L,QtXtY,guesses(i),n);
		opt.pred = bigpred_primal(Xva,yva,opt);
		opt.perf = opt.hoperf(Xva,yva,opt);
		%p{i} = perf(scores,yho,{'precrec'});
		for t = 1:T
			ap(i,t) = opt.perf.forho(t);
		end	
	end	
	[dummy,idx] = max(ap,[],1);	
	vout.lambdas = guesses(idx);
	vout.forho = ap;
	vout.guesses = guesses;
