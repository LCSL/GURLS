function [opt] = gurls_train(X,y)

%	gurls_train(X,y)
%	Trains a linear model using the dual formulation of RLS
%	The regularization parameter is chosen using a leave-one-out crossvalidation.
%	
%	INPUTS:
%		- X : 	input data matrix (one sample per row)
%		- y : 	labels (must have as many entries as there are rows in X).
%			must contain integers [1...T] with T being the number of classes.
%	OUTPUTS:
%		- opt :	struct of options with the following fields:
%			* several standard opt set by defopt.
%			* kernel 	: output of kernel_linear.m
%			* paramsel 	: output of paramsel_loocvdual.m
%			* rls		: output of rls_dual.m
%
%	The opt structure will be save to disk in a file called "quickanddirty.mat"

	T = max(y);
	codes = 2*eye(T) - 1;
	y = codes(y,:);

	name = 'quickanddirty';
	opt = defopt(name);
	opt.seq = {'kernel:linear', 'paramsel:loocvdual', 'rls:dual', 'pred:dual', 'perf:macroavg'};
	opt.process{1} = [2,2,2,0,0];
	opt.process{2} = [3,3,3,2,2];
	opt = gurls (X, y, opt,1);
