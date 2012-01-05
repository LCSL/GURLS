function [kernel] = kernel_rbf(X,y,opt)

% 	kernel_rbf(X,y,opt)
%	Computes the kernel matrix for a Gaussian kernel.
%	opt.paramsel.sigma is needed, It is set automatically by
% 	paramsel_sig* routines.

if ~isfield(opt,'kernel')
	opt.kernel.type = 'rbf';
end
if ~isfield(opt.kernel,'distance')
	opt.kernel.distance = distance(X',X');
	kernel.distance = opt.kernel.distance;
end	

%D = -(opt.kernel.distance.^2);
D = -(opt.kernel.distance);
K = exp(D/(opt.paramsel.sigma^2));

kernel.type = 'rbf';
kernel.K = K;
