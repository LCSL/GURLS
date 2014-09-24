function [kernel] = kernel_chisquared(X,y, opt)
% kernel_chisquared(opt)
% Computes the Kernel matrix for chi-squared kernel.		
%
% INPUTS:
% -OPT: must contain the field X
% 
% OUTPUT: struct with the following fields:
% -type: 'chisquared'
% -K: kernel matrix

kernel = opt.kernel;
kernel.type = 'chisquared';

if ~isfield(kernel, 'kerrange')
    kernel.kerrange = 1;
end

if ~isfield(kernel, 'init') || ~kernel.init
    for i = 1:size(X,1)
        for j = 1:i
            kernel.K(i,j) = sum(...
                            ( (X(i,:) - X(j,:)).^2 ) ./ ...
                            ( 0.5*(X(i,:) + X(j,:)) + eps));
            kernel.K(j,i) = kernel.K(i,j);
        end
    end
else
    kernel.init = 0;
end

		
