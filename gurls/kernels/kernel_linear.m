function [kernel] = kernel_linear(X, y, opt)
%	kernel_linear(X,y,opt)
%	Computes the Kernel matrix for a linear model.		

kernel.type = 'linear';
kernel.K = X*X';

