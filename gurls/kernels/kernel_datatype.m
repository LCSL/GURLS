function [kernel] = kernel_datatype(X,y, opt)
% kernel_linear(OPT)
% Load the kernel supplied in the variable X.		
%
% INPUTS:
% -OPT: not used
% 
% OUTPUT: struct with the following fields:
% -type: 'datatype'
% -K: kernel matrix

kernel = opt.kernel;
kernel.type = 'datatype';
kernel.kerrange = 1;

if ~isfield(kernel, 'init') || ~kernel.init
    kernel.K = X;
else
    kernel.init = 0;
end





