function [kernel] = kernel_Example(X, Y, opt)

% 	kernel_Example(opt)
%	Computes the kernel matrix for a Gaussian kernel.
%	INPUTS:
%		-OPT: struct with the following options:
%			- paramsel : struct containing the following fields (computed by paramsel_*).
%				- sigma : width of the gaussian kernel.
%           -X: input data matrix
%	
%	OUTPUT: struct with the following fields:
%		-type: 'rbf'
%		-K: kernel matrix

if ~isprop(opt,'kernel')
    opt.newprop('kernel', struct());
end
kernel = opt.kernel;

if ~isfield(kernel,'distance')
    if ~isfield(kernel, 'Y')
        kernel.distance = square_distance(X',X');
    else
        kernel.distance = square_distance(X',kernel.Y');
    end
end	

if ~isfield(kernel, 'kerrange')
    kernel.kerrange = {1, 2, 3}; %example
end

if ~isfield(kernel, 'init') || ~kernel.init
    
    if ~isfield(opt.paramsel, 'sigma')
        sigma = kernel.kerrange{opt.paramsel.sigmanum};
    else
        sigma = opt.paramsel.sigma{1};
    end
    
    
    %% Computation of the kernel
    D = -(kernel.distance);
    K = exp(D/(sigma^2));

    kernel.type = 'Example';
    kernel.K = K;
else
    kernel.init = 0;
end
