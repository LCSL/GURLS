function [acc] = test_classifier (W, opt)
	opt.rls.W = W;
	opt.pred = pred_primal(opt.Xte, opt.yte, opt);
	opt.perf   = perf_macroavg(opt.Xte, opt.yte, opt);
	acc = mean([opt.perf.acc]);
end
