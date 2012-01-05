function [cfr] = bigrls_primal(bX, bY, opt)

%	if (opt.kernel.primal_already_computed)
%		XtX = opt.kernel.K;
%	else
		XtX = ba_product_uvt(bX, bX);
%	end
	
	Xty = ba_product_uvt(bX, bY);
	
	n = bX.NumItems();
	
	lambda = opt.singlelambda(opt.paramsel.lambdas);
	
	cfr.W = rls_primal_driver( XtX, Xty, n, lambda);
	cfr.C = [];
	cfr.X = [];
	
end
