function [out] = conf_boltzman(X,y,opt)
% conf_boltzman(X,Y,OPT)	
% Computes the probability of belonging to the highest scoring class.
% The scores are converted in probabilities using the Boltzman distribution.
%
% INPUTS:
% -X: not used
% -Y: not used
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- pred (set by the pred_* routine)
% 
% OUTPUTS: structure with the following fields:
% - confidence: nx1 array containing  the confidence for each row of opt.pred.
% - labels: nx1 array containing predicted class for each row of opt.pred.

out = struct;
[n,k] = size(opt.pred);
expscores = exp(opt.pred);
expscores = expscores./(sum(expscores,2)*ones(1,k));
[out.confidence, out.labels] = max(expscores,[],2);
