function [cfr] = bigrls_dprimal(bX, bY, opt)

%	bigrls_dprimal(X,y,opt)
%	computes a classifier for the primal formulation of RLS.
%	The regularization parameter is set to the one found in opt.paramsel (set by the bigparamsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%	This particular implementation assumes XtX and Xty have already been computed using gdm.
%
%	NEEDS:
%		- opt.singlelamda
%		- opt.paramsel.lambdas
%		- opt.files.XtX_filename
%		- opt.files.Xty_filename

	t = load(opt.files.XtX_filename);	XtX = t.data;
	t = load(opt.files.Xty_filename);	Xty = t.data;

	n = bX.NumItems();
	
	lambda = opt.singlelambda(opt.paramsel.lambdas);
	
	cfr.W = rls_primal_driver( XtX, Xty, n, lambda);
	cfr.C = [];
	cfr.X = [];
	
end
