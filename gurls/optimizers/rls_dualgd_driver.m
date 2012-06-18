function gd = rls_dualgd_driver(X,y,opt)

% utility function called by rls_dualgd
% computes a single step for gradient descent algorithm
% NEEDS:
%   - opt.gd.method (0 for standard gd, 1 for accelerated)
%	- opt.gd.alpha1
%	- opt.gd.nu
%	- opt.gd.iter
%   - opt.paramsel.eta
%	- opt.kernel.K
%	- opt.kernel.type


gd = opt.gd;
if (opt.gd.method == 1)
	alpha1 = opt.gd.alpha1;
	nu = opt.gd.nu;
end
i = opt.gd.iter;
c = opt.gd.c;
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
gd.c = c;
if (opt.gd.method == 1)
	gd.alpha1 = alpha1;
end
