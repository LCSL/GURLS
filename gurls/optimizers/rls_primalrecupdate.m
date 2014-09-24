function [rls] = rls_primalrecupdate (X,y, opt)

% rls_primalrecupdate(X,y,opt)
% computes a classifier for the primal formulation of RLS, using a
% recursive update, starting from an initial estimator found in OPT.RLS.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- rls.W (set by the rls_primalrecinit)
%       - rls.Cinv (set by rls_primalrecinit)
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
Cinv = opt.rls.Cinv;
for i = 1:size(X,1);
    Cx = Cinv*X(i,:)';
    xCx = X(i,:)*Cx;
    Cinv = Cinv - Cx*Cx'./(1+xCx);
    W = W +(Cx*(y(i,:)-X(i,:)*W))./(1+xCx);
end

rls.W = W;
rls.C = [];
rls.X = [];
rls.Cinv = Cinv;

