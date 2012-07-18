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

sQ = size(Q,1);		
L = L + n*lambda;
L = L.^(-1);
% L = spdiags(L,0,sQ,sQ);
sL = length(L);
L = spdiags(L,0,sL,sQ);

C = (Q*L)*QtY;
