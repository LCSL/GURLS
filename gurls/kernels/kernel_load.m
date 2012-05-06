function [kernel] = kernel_load(X, y, opt)

% 	kernel_load(X,y,opt)
%	Loads the kernel matrix from disk.
%	INPUTS:
%		-X: input data matrix
%		-y: not used 
%		-OPT: struct with the following options:
%			- trainkernel : name of the file where the kernel matrix is saved.
%	
%	OUTPUT: struct with the following fields:
%		-type: 'load'
%		-K: kernel matrix
		
	if ~isfield(opt,'kernel')
		opt.kernel.type = 'load';
	end
	load(opt.trainkernel);
	kernel.type = 'load';
	kernel.K = K;
