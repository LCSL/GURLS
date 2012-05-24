function [vout] = paramsel_hoprimalgd(X, y, opt)

%	paramsel_hoprimalgd(X,y,opt)
%	Performs parameter selection when the primal formulation of gradient
% descent is used
%	The performance measure specified by opt.hoperf is maximized.
%
%	NEEDS:	
%   - opt.gd.method (0 for standard gd, 1 for accelerated)
%   - opt.gd.maxiter
%   - opt.gd.eta_numerator
%		- opt.split
%		- opt.nholdouts
%		- opt.hoperf

if (isfield(opt, 'paramsel'))
  vout = opt.paramsel;
end

[n,d] = size(X);
T = size(y,2);

% Force opt.split to be a cell array, even if just using single holdout
if ~iscell(opt.split)
  s = opt.split;
  opt = rmfield(opt, 'split');
  opt.split{1} = s;
end

% Get step size
max_eig = eigs(double(X'*X),1);
eta = opt.gd.eta_numerator/max_eig;
vout.eta = eta;
if (opt.gd.method == 1)
  nu = opt.gd.nu;
end

vout.ho_iter = zeros(1,opt.nholdouts);
for j=1:opt.nholdouts
  tr = opt.split{j}.tr;
  va = opt.split{j}.va;
  
  XtX = X(tr,:)'*X(tr,:); % d x d matrix.
  Xty = X(tr,:)'*y(tr,:); % d x T matrix.

  % Initialize
  W = zeros(d,T);
  if (opt.gd.method == 1)
    alpha1 = zeros(d,T);
  end
  
  ap = zeros(opt.gd.maxiter,T);

  % Iterate
  for i = 1:opt.gd.maxiter
    % Update coefficients
    if (opt.gd.method == 0)
      W = W + eta*(Xty - XtX*W);
    elseif (opt.gd.method == 1)
      u=((i-1)*(2*i-3)*(2*i+2*nu-1))/((i+2*nu-1)*(2*i+4*nu-1)*(2*i+2*nu-3));
      w=4*(((2*i+2*nu-1)*(i+nu-1)) /((i+2*nu-1)*(2*i+4*nu-1)) );
      alpha2 = alpha1;
      alpha1 = W;
      W = alpha1 + u*(alpha1 - alpha2) +(w*eta)*(Xty - XtX*alpha1);
    else
      error('invalid opt.gd.method')
    end
    
    % Make predictions
    opt.pred = X(va,:)*W;
    
    % Get accuracy
    opt.perf = opt.hoperf([],y(va,:),opt);
    ap(i,:) = opt.perf.acc;
    
     % To do: Add early termination condition
    
  end
  
  [~, best_idx] = max(mean(ap,2));
  vout.ho_iter(j) = best_idx;
  vout.forho{j} =ap;
  
  
end

% Average over holdouts
vout.iter = mean(vout.ho_iter);