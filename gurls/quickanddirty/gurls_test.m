function [yhat, acc] = gurls_test(X,y,opt)

%	gurls_test(X,y,opt)
%	Predicts labels on a test-set using a linear model.
%	The regularization parameter is chosen using a leave-one-out crossvalidation.
%	
%	INPUTS:
%		- X : 	input data matrix (one sample per row)
%		- y : 	labels (must have as many entries as there are rows in X).
%			must contain integers [1...T] with T being the number of classes.
%		- opt : Structure containing the fields outlined in gurls_train.m
%
%	OUTPUTS:
%		- yhat :	predicted labels (will have as many entries as there are rows in X.)
%		- acc  : 	Accuracy per class ( i.e. tp/(tp+fn)  per class ).
%		- The following fields are appended to opt:
%			* pred	: output of pred_dual.m
%			* perf	: output of perf_macroavg.m
%
%	The opt structure will be save to disk in a file called "quickanddirty.mat"

	T = max(y);
	codes = 2*eye(T) - 1;
	y = codes(y,:);

	opt = gurls (X, y, opt,2);
	[ans,yhat] = max(opt.pred,[],2);
	acc = opt.perf.acc;
