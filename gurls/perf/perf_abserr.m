function [p] = perf_abserr(X, y, opt)

% perf_rmse(X, y, opy)		
% Computes the absolute mean error for the predictions.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: structure of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%       -pred (set by the pred_* routines)
% 
% OUTPUT: struct with the following fields:
% -abserr: array of abserr for each class/output
% -forplot: ""
% -forho: array of -abserr for each class/output


if isstruct(opt.pred)
	opt.pred = opt.pred.means;
end	
if isfield (opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end
diff 		= opt.pred - y;
p.abserr		= mean(abs(diff));
p.forho 	= -p.abserr;
p.forplot 	= p.abserr;
