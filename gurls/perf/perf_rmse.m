function [p] = perf_rmse(X, y, opt)

%	 perf_rmse(X, y, opy)
%		
% 	Computes the root mean squared error for the predictions.
%
%	NEEDS:
%		opt.pred

if isfield (opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end

T 		= size(y,2);
diff 	= opt.pred - y;
p.rmse 	= norm(diff,'fro') / sqrt(numel(diff));
