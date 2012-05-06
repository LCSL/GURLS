function [out] = conf_boltzmangap(X,y,opt)

%	conf_boltzmangap(X,y,opt)	
%	computes a confidence estimation for the predicted class (i.e. highest scoring class).
%	The scores are converted in probabilities using the
%	Boltzman distribution and the difference between the highest
%	scoring class and the second highest scoring class is used as
% 	an estimate.
%
%	 INPUTS:
% 		-X: input data matrix
% 		-y: labels matrix
% 		-OPT: struct of options with the following fields:
%			- pred : predicted values (computed by the pred_* routines).
% 	OUTPUT: 
%		struct with only on field (confidence) which contains
%		one entry for each row of opt.pred.

		out.confidence = opt.pred;
		[n,k] = size(opt.pred);
		out.confidence = exp(out.confidence);
		out.confidence = out.confidence./(sum(out.confidence,2)*ones(1,k));
		out.confidence = sort(out.confidence,2,'descend');
		out.confidence = out.confidence(:,1) - out.confidence(:,2);
		
