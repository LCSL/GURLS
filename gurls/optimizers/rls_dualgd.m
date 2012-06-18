function [cfr] = rls_dualgd(X, y, opt)

%	rls_dualgd(X,y,opt)
%	computes a classifier for the dual formulation of RLS using gradient
%	descent
%	The stopping parameter is set to the one found in opt.paramsel (set by the paramsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singleiter function.
%
%	NEEDS:
%   - opt.gd.method (0 for standard gd, 1 for accelerated)
%		- opt.gd.singleiter
%		- opt.paramsel.iter
%   	- opt.paramsel.eta
%		- opt.kernel.K
%		- opt.kernel.type

iter = round(opt.gd.singleiter(opt.paramsel.iter));

%fprintf('\tSolving dual by gradient descent...\n');

n = size(opt.kernel.K,1);
T = size(y,2);

% Initialize
opt.gd.c = zeros(n,T);
if (opt.gd.method == 1)
  opt.gd.alpha1 = zeros(n,T);
end
% Iterate
for i=1:iter
	opt.gd.iter = iter;
	opt.gd = rls_dualgd_driver(X,y,opt);
end

cfr.C = opt.gd.c;
cfr.X = X;

if strcmp(opt.kernel.type, 'linear')
	cfr.W = X'*cfr.C;
	cfr.C = [];
	cfr.X = [];
else
	cfr.W = [];
end
