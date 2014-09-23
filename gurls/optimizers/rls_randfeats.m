function [rls] = rls_randfeats (X,y, opt)

% rls_primalRF(X, y, opt)
% computes a classifier for the primal formulation of RLS using the Random
% Features approach proposed in: 
%   Ali Rahimi, Ben Recht;
%   Random Features for Large-Scale Kernel Machines;
%   in Neural Information Processing Systems (NIPS) 2007.
%
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:

% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
%       - randfeat.D
%       - randfeat.samplesize
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -C: empty matrix
% -X: empty matrix

lambda = opt.singlelambda(opt.paramsel.lambdas);

%fprintf('\tSolving primal RLS...\n');

n = size(X,1);

if or(opt.randfeats.samplesize < 0, opt.randfeats.samplesize > n)
    ni = n;
else 
    ni = opt.randfeats.samplesize;
end

[XtX,Xty,rls.proj] = rp_factorize_large_real(X',y',opt.randfeats.D,'gaussian',ni);


rls.W = rls_primal_driver( XtX, Xty, n, lambda );
rls.C = [];
rls.X = [];

