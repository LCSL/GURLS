function [p] = perf_gpregr(X, y, opt)

%	perf_gpregr(X, y, opy)
% 	Computes the predictive log probability for gaussian process regression
%
% INPUTS:
% -X: not used
% -y: labels matrix
% -OPT: structure of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%       -pred (set by the pred_* routines)
% 
% OUTPUT: struct with the following field:
% -logprob: predictive log probability

if isfield (opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end

[n,T] = size(y);
diff = (opt.pred.means - y)./repmat(opt.pred.vars,1,T);
p.logprob = -(log(sum(opt.pred.vars.^2)) + norm(diff,'fro')^2 + n*log(2*pi))/2;