function [kernel] = kernel_load(X, y, opt)

% 	kernel_load(X,y,opt)
%	Loads the kernel matrix from disk.
%
%	NEEDS:
%		opt.trainkernel
		
	if ~isfield(opt,'kernel')
		opt.kernel.type = 'load';
	end
	load(opt.trainkernel);
	kernel.type = 'load';
	kernel.K = K;
