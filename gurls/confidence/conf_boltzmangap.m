function [out] = conf_boltzmangap(X,y,opt)

%	conf_boltzmangap(X,y,opt)	
%	computes the probability of belonging to the highest scoring class.
% 	Output:
%		- out: 
%		struct with only on field (confidence) which contains
%		one entry for each row of opt.pred.
%
%	The scores are converted in probabilities according to the
%	Boltzman distribution and the difference between the highest
%	scoring class and the second highest scoring class is used as
% 	an estimate.
%
% NEEDS:
%		- opt.pred


		out.confidence = opt.pred;
		[n,k] = size(opt.pred);
		out.confidence = exp(out.confidence);
		out.confidence = out.confidence./(sum(out.confidence,2)*ones(1,k));
		out.confidence = sort(out.confidence,2,'descend');
		out.confidence = out.confidence(:,1) - out.confidence(:,2);
		
