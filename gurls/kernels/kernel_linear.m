function [kernel] = kernel_linear(X,y, opt)
% kernel_linear(OPT)
% Computes the Kernel matrix for a linear model.		
%
% INPUTS:
% -OPT: must contain the field X
% 
% OUTPUT: struct with the following fields:
% -type: 'linear'
% -K: kernel matrix
kernel = opt.kernel;
kernel.type = 'linear';
kernel.kerrange = 1;

if ~isfield(opt.kernel, 'init') || ~opt.kernel.init
    kernel.K = X*X';
else
    kernel.init = 0;
end


