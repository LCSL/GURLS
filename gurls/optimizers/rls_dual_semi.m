function [cfr] = rls_dual_semi (X,y, opt)
% rls_dualsemi(X, y, opt)
% computes a classifier for semi-supervised learning with the dual 
% formulation of RLS.
% The parameters for both the manifold and ambient regularization terms 
% are set in opt.paramsel.
%
% INPUTS:
% -OPT: struct of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_siglam_semi routine)
%		- paramsel.lambdas_m (set by the paramsel_siglam_semi routine)
%		- paramsel.sigma_m (set by the paramsel_siglam_semi routine)
%       - kernel (set by the kernel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
%		- kernel.type
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% if kernel.type='linear'
% -W: matrix of coefficient vectors of primal rls estimator for each class
% -C: empty matrix
% -X: empty matrix
% else
% -W: empty matrix
% -C: matrix of coefficient vectors of dual rls estimator for each class
% -X: empty matrix

lambda = opt.singlelambda(opt.paramsel.lambdas);
lambda_m = opt.singlelambda(opt.paramsel.lambdas_m);
sigma_m = opt.paramsel.sigma_m;


n = size(opt.kernel.K,1);
T = size(y,2);

% find range of labelled y
labelled = (~isnan(y(:,1)));
J = diag(+labelled);
m = sum(+labelled);

% find graph laplacian G
D = -(opt.kernel.distance);
W = exp(D/(sigma_m^2)).*(opt.kernel.distance < 4.*sigma_m^2);
G = diag(sum(W,2))-W;

Jy = zeros(size(y));
Jy(labelled,:) = y(labelled,:);

try
    cfr.C =(J*opt.kernel.K + (m*lambda)*eye(n) + (m*lambda_m)*G*opt.kernel.K) \ Jy;
	%R = chol(J*opt.kernel.K + (m*lambda)*eye(n) + (m*lambda_m)*G*opt.kernel.K); 
	%cfr.C = R\(R'\Jy);
catch
    warning('using SVD');
	[Q,L,V] = svd(J*opt.kernel.K + (m*lambda_m)*G*opt.kernel.K);
	Q = double(Q);
	L = double(diag(L));
	cfr.C = rls_eigen(Q,L,Q'*Jy,lambda*m/n,n);
end	

if strcmp(opt.kernel.type, 'linear')
	cfr.W = X'*cfr.C;
	cfr.C = [];
	cfr.X = [];
else
	cfr.W = [];
	cfr.X = X;
end

end

