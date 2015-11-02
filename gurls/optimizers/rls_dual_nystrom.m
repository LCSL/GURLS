function [cfr] = rls_dual_nystrom(X,y, opt)

% rls_dual_nystrom(X, y, opt)
% computes a predictor for Nystrom KRLS.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -OPT: struct of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routine)
%       - kernel.Knm (set by the kernel_nystrom routines)
%       - kernel.Kmm (set by the kernel_nystrom routines)
%   fields with default values set through the defopt function:
%		- singlelambda
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: empty matrix
% -C: matrix of coefficient vectors of dual rls estimator for each class
% -X: training samples used by the routine

lambda = opt.singlelambda(opt.paramsel.lambdas);

n = size(opt.kernel.Knm,1);

cfr.C = ( opt.kernel.Knm' * opt.kernel.Knm + n * lambda * opt.kernel.Kmm) \ opt.kernel.Knm' * y;

% TO DO: support for linear kernel
% if strcmp(opt.kernel.type, 'linear')
% 	cfr.W = X'*cfr.C;
% 	cfr.C = [];
% 	cfr.X = [];
% else
	cfr.W = [];
	cfr.X = X;
end
