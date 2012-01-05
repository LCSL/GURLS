function [out] = conf_maxscore(X,y,opt)

%	conf_maxscore(X,y,opt)	
%	computes the probability of belonging to the highest scoring class.
% 	Output:
%		- out: 
%		struct with only on field (confidence) which contains
%		one entry for each row of opt.pred.
%
%	NEEDS:
%		-opt.pred

		out = struct;
		[out.confidence, out.labels] = max(opt.pred,[],2);

