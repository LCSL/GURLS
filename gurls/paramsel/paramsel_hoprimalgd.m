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
%	- opt.split
%	- opt.nholdouts
%	- opt.hoperf

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
opt.paramsel.eta = eta;
if (opt.gd.method == 1)
  	nu = opt.gd.nu;
end

vout.ho_iter = zeros(1,opt.nholdouts);
for j=1:opt.nholdouts
  	tr = opt.split{j}.tr;
  	va = opt.split{j}.va;

	opt.gd.W = zeros(d,T);
  	if (opt.gd.method == 1)
    	opt.gd.alpha1 = zeros(d,T);
  	end
  	ap = zeros(opt.gd.maxiter,T);

  	% Iterate
  	for i = 1:opt.gd.maxiter
		opt.gd.iter = i;
		opt.gd = rls_primalgd_driver(X(tr,:),y(tr,:),opt);
		opt.rls.W = opt.gd.W;
    	% Make predictions
    	opt.pred = pred_primal(X(va,:),[],opt);

    
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
