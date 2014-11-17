function [rls] = rls_primalrecupdatecholesky (X, y, opt)

% rls_primalrecupdatecholesky(X,y,opt)
% computes a classifier for the primal formulation of RLS, using a
% recursive update of the right Cholesky factor R of covariance matrix C = XtX,
% starting from an initial R found in OPT.RLS.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: Options object with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- rls.W (set by rls_primalrecinitcholesky)
%       - rls.R (set by rls_primalrecinitcholesky)
%       - rls.b (set by rls_primalrecinitcholesky)
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -Cinv: inverse of the regularized kernel matrix in the primal space
% -C: empty matrix
% -X: empty matrix

W = opt.rls.W;
R = opt.rls.R;
b = opt.rls.b;

n = size(X,1);

% Sequence of rank-1 updates by application of Cholesky rank-1 updates
for i = 1:n;
    
    % Update b
    b = b + X(i,:)'*y(i,:);
    
    % Update Cholesky factor R
    R = cholupdate(R,X(i,:)');

end

W = R\(R'\b);

rls.b = b;
rls.W = W;
rls.C = [];
rls.X = [];
rls.XtX = [];
rls.R = R;

