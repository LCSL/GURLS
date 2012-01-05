function [cfr] = rls_monline_driver(X,y,opt)


lambda = opt.singlelambda(opt.paramsel.lambdas);

fprintf('\tSolving primal RLS...\n');

[n,d] = size(X);

XtX = X'*X; % d x d matrix.
Xty = X'*y + n*lambda*opt.rls.W0; % d x T matrix.

cfr.W = rls_primal_driver( XtX, Xty, n, lambda );
cfr.C = [];
cfr.X = [];

