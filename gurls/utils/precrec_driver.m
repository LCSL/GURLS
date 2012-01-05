function [ap] = precrec_driver(out, gt, draw)

% compute precision/recall
[so,si]=sort(-out);
tp=gt(si)>0;
fp=gt(si)<0;

fp=cumsum(fp);
tp=cumsum(tp);
rec=tp/sum(gt>0);
prec=tp./(fp+tp);
% compute average precision

ap=0;
for t=0:0.1:1
    p=max(prec(rec>=t));
    if isempty(p)
        p=0;
    end
    ap=ap+p/11;
end

if draw
    % plot precision/recall
	figure();
    plot(rec,prec,'-');
    grid;
    xlabel 'recall'
    ylabel 'precision'
    %title(sprintf('class: %s, subset: %s, AP = %.3f',cls,VOCopts.testset,ap));
end
