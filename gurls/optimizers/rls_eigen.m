function C = rls_eigen(Q,L,QtY,lambda,n)

sQ = size(Q,1);		
L = L + n*lambda;
L = L.^(-1);
L = spdiags(L,0,sQ,sQ);

%C = Q*L*Q'*Y;
C = (Q*L)*QtY;
