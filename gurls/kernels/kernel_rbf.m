function [kernel] = kernel_rbf(X,y,opt)

% 	kernel_rbf(X,y,opt)
%	Computes the kernel matrix for a Gaussian kernel.
%	INPUTS:
%		-X: input data matrix
%		-y: not used 
%		-OPT: struct with the following options:
%			- paramsel : struct containing the following fields (computed by paramsel_*).
%				- sigma : width of the gaussian kernel.
%	
%	OUTPUT: struct with the following fields:
%		-type: 'rbf'
%		-K: kernel matrix


if ~isfield(opt,'kernel')
	opt.kernel.type = 'rbf';
end


if ~isfield(opt.kernel,'distance')
	opt.kernel.distance = distance(X',X');
end	

kernel = opt.kernel;
D = -(opt.kernel.distance);
K = exp(D/(opt.paramsel.sigma^2));

kernel.type = 'rbf';
kernel.K = K;
