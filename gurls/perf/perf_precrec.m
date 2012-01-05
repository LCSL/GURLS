function [p] =  perf_precrec(X, y, opt)

%	perf_precrec(X,y,opt)
% 	Computes the average precision per class
%	NEEDS:
%		- opt.pred
%


% Code adapted from PASCAL VOC 2007 code

if isfield (opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end

y_true = y;
y_pred = opt.pred;

T = size(y,2);
for t = 1:T,
	p.ap(t) = precrec_driver(y_pred(:,t), y_true(:,t),0);
	p.forho(t) = p.ap(t);
	p.forplot(t) = p.ap(t);
end

