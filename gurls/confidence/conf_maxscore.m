function [out] = conf_maxscore(X,y,opt)

%	conf_maxscore(X,y,opt)	
%	computes a codfidence estimation for the predicted class (i.e. highest scoring class).
%	The difference between the highest scoring class and the second highest
%	scoring class is considered.
%
%	INPUT:
%		- X 	: input data matrix (not used).
%		- y 	: labels matrix (not used).
%		- opt 	: struct of options containing the following fields:
%				- pred : predicted valies (computed by the pred_* routines).
% 	Output:
%		struct with fields:
%			- confidence 	: confidence estimate for the predicted label.
%			- labels	: predicted label.

		out = struct;
		[out.confidence, out.labels] = max(opt.pred,[],2);

