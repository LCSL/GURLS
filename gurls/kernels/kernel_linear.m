function [kernel] = kernel_linear(X, y, opt)
% kernel_linear(X,Y,OPT)
% Computes the Kernel matrix for a linear model.		
%
% INPUTS:
% -X: input data matrix
% -Y: not used 
% -OPT: not used
% 
% OUTPUT: struct with the following fields:
% -type: 'linear'
% -K: kernel matrix

kernel.type = 'linear';
kernel.K = X*X';

