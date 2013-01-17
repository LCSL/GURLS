function [cfr] = bigrls_primal(bX, bY, opt)

%	bigrls_primal(X,y,opt)
%	computes a classifier for the primal formulation of RLS.
%	The regularization parameter is set to the one found in opt.paramsel (set by the bigparamsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.

%	INPUT:
%		- X : input data bigarray
%		- Y : labels bigarray
%		- OPT : struct witht he following fields:
%			- Fields set by other biggruls tasks:
%				* paramsel.lambdas (set by the bigparamsel_*) routines.
%			- Fields set through the bigdefopt function:
%				* singlelambda
%
%	OUTPUT: structure with the following fields:
%		- W : matrix of coefficient vectors of rls estimator for each class
%		- C : empty matrix
%		- X : empty matrix

	
	Xty = ba_product_uvt(bX, bY);
	
	n = bX.NumItems();
	
	lambda = opt.singlelambda(opt.paramsel.lambdas);
	
	cfr.W = rls_primal_driver( XtX, Xty, n, lambda);
	cfr.C = [];
	cfr.X = [];
	
end
