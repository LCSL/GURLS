function [vout] = bigparamsel_hoprimal(X,y,opt)

%	bigparamsel_hoprimal(X,y,opt)
%	Performs parameter selection when the primal formulation of RLS is used 
%	This function uses a hold-out set to perform parameter selection.
%	The performance measure specified by opt.hoperf is maximized.
%
%	INPUT:
%		- X : input data bigarray
%		- Y : labels bigarray
%		- opt : struct witht he following fields:
%			- Fields set through the bigdefopt function:
%				* nlambdaa
%				* hoperf
%			- Fields that need to be set by hand:
%
%				* opt.files.Xva_filename 	: Validation data bigarray
%				* opt.files.yva_filename 	: Validation labels bigarray
%
%	OUTPUT: structure with the following fields:
%		- lambdas 	: values of the regulariazation parameter maximizing the
%			    	  validation performance (one for each class).	
%		- forho		: matrix with validation performance for each class and for each value of the paramter.
%		- guesses	: contains the values tried for the lambda paramter.



	X.Transpose(true);
	y.Transpose(true);



	Xva = bigarray.Obj(opt.files.Xva_filename);
	yva = bigarray.Obj(opt.files.yva_filename);

	Xva.Transpose(true);
	yva.Transpose(true);

	XtX = ba_product_uvt(X, X); 
	Xty = ba_product_uvt(X, y);
	XvatXva = ba_product_uvt(Xva, Xva);
	Xvatyva = ba_product_uvt(Xva, yva);


	
	
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
