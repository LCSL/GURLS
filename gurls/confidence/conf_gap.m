function [out] = conf_gap(X,y,opt)

%	conf_gap(X,y,opt)	
%	computes the probability of belonging to the highest scoring class.
% 	Output:
%		- out: 
%		struct with only on field (confidence) which contains
%		one entry for each row of opt.pred.
%	
%	The scores are converted in probabilities and then the
%	difference between the highest scoring class and the second highest
%	scoring class is considered.
%
%	NEEDS:
%		- opt.pred

		out.confidence = opt.pred;
		out.confidence = sort(out.confidence,2,'descend');
		out.confidence = out.confidence(:,1) - out.confidence(:,2);
