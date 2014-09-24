function [p] = perf_rmse(X,y, opt)

% perf_rmse(opt)		
% Computes the root mean squared error for the predictions.
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%   -y: labels matrix
%   fields that need to be set through previous gurls tasks:
%       -pred (set by the pred_* routines)
% 
% OUTPUT: struct with the following fields:
% -rmse: array of rmse for each class/output
% -forplot: ""
% -forho: array of -rmse for each class/output

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
p.rmse		= sqrt(sum(diff.^2,1)/size(diff,1));
p.forho 	= -p.rmse;
p.forplot 	= p.rmse;
