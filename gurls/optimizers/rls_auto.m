function [cfr] = rls_auto(X, y, opt)
% rls_auto(X,Y,OPT)
% computes a RLS classifier, with automatic selection of primal/dual procedure.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: structure of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routine)
%       - kernel (set by the kernel_* routines)
%		- kernel.type (needs to be set to 'linear')
%   fields with default values set through the defopt function:
%		- singlelambda
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: structure with the following fields:
% -W: matrix of coefficient vectors of primal rls estimator for each class
% -C: empty matrix
% -X: empty matrix

[n,d] = size(X);

if (n > d) % Do primal
	cfr = rls_primal(X, y, opt);
else % Do dual
	cfr = rls_dual(X, y, opt);
end
	
end
