function [cfr] = rls_dual_mkl (X,y, opt)

% rls_dual_mkl(X, y, opt)
% computes a predictor based on MKL framework with elastic net regularization.
% 
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -OPT: struct of options with the following fields (and subfields):
%   fields that need to be set manually:
%		- mkl.kernel: a cell array of pre-computed kernels
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: empty matrix
% -C: matrix of coefficient vectors of dual rls estimator for each class
% -X: training samples used by the routine

lambda = opt.singlelambda(opt.paramsel.lambdas);


n = size(opt.kernel.K,1);
T = size(y,2);

%fprintf('\tSolving dual RLS...(n = %d, % = %d)', n, T);

try
	R = chol(opt.kernel.K + (n*lambda)*eye(n)); 
	
	cfr.C = R\(R'\y);
catch
	[Q,L,V] = svd(opt.kernel.K);
	Q = double(Q);
	L = double(diag(L));
	cfr.C = rls_eigen(Q,L,Q'*y,lambda,n);
end	

if strcmp(opt.kernel.type, 'linear')
	cfr.W = X'*cfr.C;
	cfr.C = [];
	cfr.X = [];
else
	cfr.W = [];
	cfr.X = X;
end
