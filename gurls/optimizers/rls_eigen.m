function C = rls_eigen(Q,L,QtY,lambda,n)
% rls_eigen(Q,L,QTY,LAMBDA,N)
% Computes RLS estimator given the singular value decomposition of the
% kernel matrix
% 
% INPUTS
% -Q: eigenvectors of the kernel matrix
% -L: eigenvalues of the kernel matrix
% -QTY: result of the matrix multiplication of the transpose of Q times the
%       labels vector Y (Q'*Y)
% -LAMBDA: regularization parameter
% -N: number of training samples
% 
% OUTPUT:
% -C: rls coefficient vector

L = L + n*lambda;
L = L.^(-1);

% sL = length(L);
% L = spdiags(L,0,sL,sL);
L = diag(L);

C = (Q*L)*QtY;
