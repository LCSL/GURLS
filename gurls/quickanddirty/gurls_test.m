function [yhat, acc] = gurls_test(X,y,opt)

	T = max(y);
	codes = 2*eye(T) - 1;
	y = codes(y,:);

	opt = gurls (X, y, opt,2);
	[~,yhat] = max(opt.pred,[],2);
	acc = opt.perf.acc;
