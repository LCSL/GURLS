function [out] = conf_gap(X, y, opt)
% conf_gap(OPT)	
% Computes the probability of belonging to the highest scoring class.
% The scores are converted in probabilities and then the
% difference between the highest scoring class and the second highest
% scoring class is considered.

%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- pred (set by the pred_* routine)
% 
% OUTPUTS: structure with only one field (confidence) which contains
%	one entry for each row of opt.pred.
	
out.confidence = opt.pred;
out.confidence = sort(out.confidence,2,'descend');
out.confidence = out.confidence(:,1) - out.confidence(:,2);
