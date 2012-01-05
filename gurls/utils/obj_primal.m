function [obj] = obj_primal(W, X, Y, lambda)

	[n,d] = size(X);

	obj1 = norm(W,'fro');
	obj1 = lambda*obj1*obj1;

	obj2 = norm(Y - X*W,'fro');
	obj2 = obj2*obj2/n;

	obj = 0.5*(obj1 + obj2);
	
end	

