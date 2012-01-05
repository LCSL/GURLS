function [cfr] = rls_tax_d(X,Y,opt)

K=opt.kernel.K;
[V,D]=eig(K);
T=size(Y,2);
n=size(X,1);
p=size(X,2);
A=1/T*eye(T);
newA=1/T*eye(T);
I= eye(n);
code=eye(T);
lambda=opt.lambda;

for iter= 1:opt.Maxiter
A=newA;    
[U,S] =eig(A);
Sigma=diag(S);
Yrotated=YEncoded*U;
Z=zeros(n,T);
Zthilde=zeros(T,n);
 %find W ,A is fixed
 for i=1:T
BBB=diag(D+I*(lambda/Sigma(i)));  
BB=1./BBB;
B=diag(BB);
Zthilde(i,:)=V*B*V'*Yrotated(:,i);
Z=Z+Zthilde(i,:)'*U(:,i)';
 end
%find A, Z is fixed
sym=Z'*K*Z;
[base,value]=eig(sym);
for i=1:T
    if ((value(i,i)<0)||(value(i,i)==0))
        value(i,i)=0;
    end
end
newA=base*sqrt(value)*base';
newA=newA./trace(newA);
end
C=Z*pinv(A);
cfr.C=C;
cfr.X=X;
cfr.A=A;
cfr.code=code;
end
