function [kernel] = kernel_load(X,y, opt)

% 	kernel_load(opt)
%	Loads the kernel matrix from disk.
%	INPUTS:
%		-OPT: struct with the following options:
%			- trainkernel : name of the file where the kernel matrix is saved.
%	
%	OUTPUT: struct with the following fields:
%		-type: 'load'
%		-K: kernel matrix
    
kernel = opt.kernel;
kernel.type = 'load';

if ~isfield(kernel, 'kerrange')
    kernel.kerrange = 1;
end

if ~isfield(kernel, 'init') || ~kernel.init
    load(opt.trainkernel);
    kernel.K = K;
else
    kernel.init = 0;
end
