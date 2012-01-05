function [out] = conf_boltzman(X,y,opt)

%	conf_boltzman(X,y,opt)	
%	computes the probability of belonging to the highest scoring class.
% 	Output:
%		- out: 
%		struct with only on field (confidence) which contains
%		one entry for each row of opt.pred.
%	
%	The scores are converted in probabilities using the Boltzman distribution.
%
%	NEEDS:
%		- opt.pred
		
		out = struct;
		[n,k] = size(opt.pred);
		expscores = exp(opt.pred);
		expscores = expscores./(sum(expscores,2)*ones(1,k));
		[out.confidence, out.labels] = max(expscores,[],2);
