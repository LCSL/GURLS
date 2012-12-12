function [kernel] = kernel_randfeats(X, y, opt)
% kernel_linear(X,Y,OPT)
% Computes the Kernel matrix for a linear model.		
% Instead of using real features, the kernel is based on random
% projections. See the following paper for details:
%   Ali Rahimi, Ben Recht;
%   Random Features for Large-Scale Kernel Machines;
%   in Neural Information Processing Systems (NIPS) 2007.
%
% INPUTS:
% -X: input data matrix
% -Y: not used 
% -OPT: not used
% 
% OUTPUT: struct with the following fields:
% -type: 'linear'
% -K: kernel matrix

kernel.type = 'randfeats';

[n,d] = size(X);

if or(opt.randfeats.samplesize < 0, opt.randfeats.samplesize > n)
    ni = n;
else 
    ni = opt.randfeats.samplesize;
end

[XtX,Xty,W] = rp_factorize_large_real(X',y',opt.randfeats.D,'gaussian',ni);


kernel.XtX = XtX;
kernel.Xty = Xty;
kernel.W = W;