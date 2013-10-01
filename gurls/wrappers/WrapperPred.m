function y = WrapperPred(X,opt)

if isfield(opt,'C')
    opt.kernel = opt.kernelfun([],[],opt);
    opt.rls.X = opt.X;
    opt.rls.C = opt.C;
    opt.predkernel = predkernel_traintest(X,[],opt);
    y = pred_dual(X,[],opt);
elseif isfield(opt,'W')
    opt.rls.W = opt.W;
    y = pred_primal(X,[],opt);
else
    error('option structure doesn not contain estimator')
end