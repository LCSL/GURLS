function [cfr] = rls_primalgd(X, y, opt)

%	rls_primalgd(X,y,opt)
%	computes a classifier for the primal formulation of RLS using gradient
%	descent
%	The stopping parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singleiter function.
%
%	NEEDS:
%   - opt.gd.method (0 for standard gd, 1 for accelerated)
%		- opt.gd.singleiter
%		- opt.paramsel.iter
%   - opt.paramsel.eta

iter = round(opt.gd.singleiter(opt.paramsel.iter));

%fprintf('\tSolving primal by gradient descent...\n');

[n,d] = size(X);
T = size(y,2);


% Initialize
opt.gd.W = zeros(d,T);
if (opt.gd.method == 1)
  opt.gd.alpha1 = zeros(d,T);
end

% Iterate
for i=1:iter
	opt.gd.iter = iter;
	opt.gd = rls_primalgd_driver(X,y,opt);
end

cfr.W = opt.gd.W;
cfr.C = [];
cfr.X = [];

