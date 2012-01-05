function [scores] = pred_dual(X, y, opt)
if strcmp(opt.kernel.type , 'linear')
	scores = pred_primal(X, y, opt);
else
	% Write the dual prediction code here.
	scores = opt.predkernel.K*opt.rls.C;
end
