function [cfr] = rls_baglcload(X,y,opt)

%% Split X into opt.nSplits sets
nSplits = opt.nSplits;
% Use each classifier to predict on the whole dataset and end up with a Nx(TxnSplits) matrix wich we can now use to train a new classifier.
%cfr = rls_bagging(X,y,opt);
t = load(opt.bagcfrfile);
cfr = t.opt.rls;
clear t;
F = [];
t_opt = opt;
for splitno = 1:nSplits
    % Make this more general
	t_opt.rls.W = cfr.W{splitno};
	F = [F pred_primal(X,y,t_opt)];
end
opt.seq = opt.baglcseq;
opt.process = opt.baglcprocess;
gurls(F,y,opt,1);
t = load([opt.savefile]);
cfr.lc.W = t.opt.rls.W;

