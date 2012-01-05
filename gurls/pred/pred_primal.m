function [scores] = pred_primal(X, y, opt)

%	pred_primal(X,y,opt)
%	computes the predictions of the linear classifier stored in opt.rls.W (generated
%	by some rls_* method) on the points passed in the X matrix.
%
%	NEEDS:
%		- opt.rls.W

	scores = X*opt.rls.W;	
