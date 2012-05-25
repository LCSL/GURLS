function [vout] = paramsel_hodualgd(X, y, opt)

%	paramsel_hodualgd(X,y,opt)
%	Performs parameter selection when the dual formulation of gradient
% descent is used
%	The performance measure specified by opt.hoperf is maximized.
%
%	NEEDS:	
%   - opt.gd.method (0 for standard gd, 1 for accelerated)
%   - opt.gd.maxiter
%		- opt.split
%		- opt.nholdouts
%		- opt.kernel.type
%		- opt.kernel.K
%		- opt.hoperf

if (isfield(opt, 'paramsel'))
	vout = opt.paramsel;
end

[n, T] = size(y); % number of samples, number of classes for multiclass

% Force opt.split to be a cell array, even if just using single holdout
if ~iscell(opt.split)
	s = opt.split;
	opt = rmfield(opt, 'split');
	opt.split{1} = s;
end

% Get step size
if (opt.gd.method == 0)
	max_eig = eigs(double(opt.kernel.K),1);
	eta = opt.gd.eta_numerator/max_eig;
elseif (opt.gd.method == 1)
	eta = 1/n;
	nu = opt.gd.nu;
else
	error('invalid opt.gd.method')
end
vout.eta = eta;

vout.ho_iter = zeros(1,opt.nholdouts);
for j=1:opt.nholdouts
	tr = opt.split{j}.tr;
	va = opt.split{j}.va;
	  
	% Initialize
	coeffs = zeros(length(tr),T);
	if (opt.gd.method == 1)
		alpha1 = zeros(length(tr),T);
	end
	
	ap = zeros(opt.gd.maxiter,T);
	
	% Iterate
	for i = 1:opt.gd.maxiter
	  % Update coefficients
	  	if (opt.gd.method == 0)
	   		coeffs = coeffs + eta*(y(tr,:) - opt.kernel.K(tr,tr)*coeffs);
	  	elseif (opt.gd.method == 1)
	    	u=((i-1)*(2*i-3)*(2*i+2*nu-1))/((i+2*nu-1)*(2*i+4*nu-1)*(2*i+2*nu-3));
	    	w=4*(((2*i+2*nu-1)*(i+nu-1)) /((i+2*nu-1)*(2*i+4*nu-1)) );
	    	alpha2 = alpha1;
	    	alpha1 = coeffs;
	    	coeffs = alpha1 + u*(alpha1 - alpha2) +(w*eta)*(y(tr,:)- opt.kernel.K(tr,tr)*alpha1);
	  	else
	    	error('invalid opt.gd.method')
	  	end
	  
	  	% Make predictions
	  	opt.pred = opt.kernel.K(va,tr)*coeffs;
	  
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
