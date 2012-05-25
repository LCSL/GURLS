function [cfr] = rls_dualcg(X, y, opt)

%	rls_dualcg(X,y,opt)
%	computes a classifier for the dual formulation of RLS using the conjugate gradient descent
%	descent
%	The stopping parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singleiter function.
%
%	NEEDS:
%   - opt.gd.method (0 for standard gd, 1 for accelerated)
%		- opt.cg.singleiter
%		- opt.paramsel.iter
%   	- opt.paramsel.eta
%		- opt.kernel.K
%		- opt.kernel.type

[n,T] = size(y);
iter = round(opt.cg.singleiter(opt.paramsel.iter));

%fprintf('\tSolving dual by gradient descent...\n');

n = size(opt.kernel.K,1);
T = size(y,2);

K = (1/n) * opt.kernel.K;
% Initialize

c = zeros(n,T);
	
opt.cg.a = zeros(n,T);
opt.cg.r = y;
opt.cg.d = y;
opt.cg.t = K*y;

for i = 1:iter
	opt.cg = rls_dualcg_driver(X,y,opt);
end

cfr.C = n*opt.cg.a;
cfr.X = X;

if strcmp(opt.kernel.type, 'linear')
	cfr.W = X'*cfr.C;
	cfr.C = [];
	cfr.X = [];
else
	cfr.W = [];
end
