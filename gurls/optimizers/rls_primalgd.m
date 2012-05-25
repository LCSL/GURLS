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

XtX = X'*X; % d x d matrix.
Xty = X'*y; % d x T matrix.

% Initialize
W = zeros(d,T);
if (opt.gd.method == 1)
  alpha1 = zeros(d,T);
  nu = opt.gd.nu;
end

% Iterate
for i=1:iter
  if (opt.gd.method == 0)
    W = W + opt.paramsel.eta*(Xty - XtX*W);
  elseif (opt.gd.method == 1)
    u=((i-1)*(2*i-3)*(2*i+2*nu-1))/((i+2*nu-1)*(2*i+4*nu-1)*(2*i+2*nu-3));
    w=4*(((2*i+2*nu-1)*(i+nu-1)) /((i+2*nu-1)*(2*i+4*nu-1)) );
    alpha2 = alpha1;
    alpha1 = W;
    W = alpha1 + u*(alpha1 - alpha2) +(w*opt.paramsel.eta)*(Xty - XtX*alpha1);
  else
    error('invalid opt.gd.method')
  end
end

cfr.W = W;
cfr.C = [];
cfr.X = [];

