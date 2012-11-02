function [p] = perf_rmse(X, y, opt)

%	 perf_rmse(X, y, opy)
%		
% 	Computes the root mean squared error for the predictions.
%
%	NEEDS:
%		opt.pred

if isstruct(opt.pred)
	opt.pred = opt.pred.means;
end	
if isfield (opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end
T 		= size(y,2);
n 		= size(y,1);
diff 		= opt.pred - y;
p.rmse		= sqrt(sum(diff.^2,1));
p.forho 	= -p.rmse;
p.forplot 	= p.rmse;
