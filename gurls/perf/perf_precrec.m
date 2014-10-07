function [p] =  perf_precrec(X,y, opt)
% perf_precrec(opt)
% Computes the average precision per class
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%   -y: labels matrix
%   fields that need to be set through previous gurls tasks:
%       -pred (set by the pred_* routines)
% 
% OUTPUT: struct with the following fields:
% -acc: array of prediction accuracy for each class
% -forho: ""
% -forplot: ""
% 
% Code adapted from PASCAL VOC 2007 code

if isstruct(opt.pred)
	opt.pred = opt.pred.means;
end	
if isprop(opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end

y_true = y;
y_pred = opt.pred;

T = size(y,2);
if T == 1
	p.ap = precrec_driver(y_pred,y_true,0);
	p.forho = p.ap;
	p.forplot = p.ap;
else	
	for t = 1:T,
		p.ap(t) = precrec_driver(y_pred(:,t), y_true(:,t),0);
		p.forho(t) = p.ap(t);
		p.forplot(t) = p.ap(t);
	end
end

