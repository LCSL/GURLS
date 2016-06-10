function [p] = perf_rmsestd(X, y, opt)

% perf_rmse(opt)		
% Computes the root mean squared error for the predictions.
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%   -X: not used
%   -y: labels matrix
%   fields that need to be set through previous gurls tasks:
%       -pred (set by the pred_* routines)
% 
% OUTPUT: struct (opt.pred) with the following fields:
% -rsquare: array of rsquare for each class/output
% -forplot: ""
% -forho: array of -rsquare for each class/output

if isstruct(opt.pred)
	opt.pred = opt.pred.means;
end	
if isprop(opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end
T           = size(y,2);
n           = size(y,1);
diff 		= opt.pred - y;
p.rmsestd	= sqrt(sum(diff.^2,1)/sum(y.^2,1));
p.forho 	= -p.rmsestd;
p.forplot 	= p.rmsestd;
