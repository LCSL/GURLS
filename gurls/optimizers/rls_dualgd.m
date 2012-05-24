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
%   - opt.paramsel.eta
%		- opt.kernel.K
%		- opt.kernel.type

iter = round(opt.gd.singleiter(opt.paramsel.iter));

fprintf('\tSolving dual by gradient descent...\n');

n = size(opt.kernel.K,1);
T = size(y,2);

% Initialize
c = zeros(n,T);
if (opt.gd.method == 1)
  alpha1 = zeros(n,T);
  nu = opt.gd.nu;
end

% Iterate
for i=1:iter
  if (opt.gd.method == 0)
    c = c + opt.paramsel.eta*(y - opt.kernel.K*c);
  elseif (opt.gd.method == 1)
    u=((i-1)*(2*i-3)*(2*i+2*nu-1))/((i+2*nu-1)*(2*i+4*nu-1)*(2*i+2*nu-3));
    w=4*(((2*i+2*nu-1)*(i+nu-1)) /((i+2*nu-1)*(2*i+4*nu-1)) );
    alpha2 = alpha1;
    alpha1 = c;
    c = alpha1 + u*(alpha1 - alpha2) +(w*opt.paramsel.eta)*(y - opt.kernel.K*alpha1);
  else
    error('invalid opt.gd.method')
  end
end

cfr.C = c;
cfr.X = X;

if strcmp(opt.kernel.type, 'linear')
	cfr.W = X'*cfr.C;
	cfr.C = [];
	cfr.X = [];
else
	cfr.W = [];
end