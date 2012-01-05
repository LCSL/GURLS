function [scores] = pred_primal(X, y, opt)

t_opt = struct;
t_opt.rls = struct;


%% 1. averaged score
F = [];
for i = 1:numel(opt.rls.W) % for each of the bagging classifier
	t_opt.rls.W = opt.rls.W{i};
	F = [F pred_primal(X,y,t_opt)];
end
scores = F*opt.rls.lc.W;
