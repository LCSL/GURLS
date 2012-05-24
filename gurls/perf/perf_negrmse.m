function [p] = perf_rmse(X, y, opt)


%	 perf_negrmse(X, y, opy)
%		
% 	Computes the negative root mean squared error for the predictions.
%	The negative is used to match the "max" operation in paramesel_* routines.
%
%	NEEDS:
%		opt.pred


p = perf_rmse(X,y,opt);
p.rmse = -p.rmse;
