function [W] = rls_primal_driver(XtX, Xty, n, lambda)
d = size(XtX,1);

try % Cholesky
	XtX = XtX + (n*lambda)*eye(d);
	R = chol(XtX);
	W = R\(R'\Xty);

catch  % SVD
	[Q,L,V] = svd(XtX);
	Q = double(Q);
	L = double(diag(L));
	QtXtY = Q'*Xty;
	
	% regularization is done inside rls_eigen
	W = rls_eigen(Q, L, QtXtY, lambda, n);
end



end
