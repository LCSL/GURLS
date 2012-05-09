function [opt] = gurls_train(X,y)

	T = max(y);
	codes = 2*eye(T) - 1;
	y = codes(y,:);

	name = 'quickanddirty';
	opt = defopt(name);
	opt.seq = {'kernel:linear', 'paramsel:loocvdual', 'rls:dual', 'pred:dual', 'perf:macroavg'};
	opt.process{1} = [2,2,2,0,0];
	opt.process{2} = [3,3,3,2,2];
	opt = gurls (X, y, opt,1);
