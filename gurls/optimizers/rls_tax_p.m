function [cfr ]  = rls_tax_p(X,Y,opt)

X = X';

[V,D] = eig(X*X');

T = size(Y,2);
n = size(X,2);
p = size(X,1);
A = 1/T*eye(T);
newA = 1/T*eye(T);

I =  eye(p);

code = eye(T);
Wthilde = zeros(T,p);

YEncoded = 0.5*(Y + 1);

lambda = opt.singlelambda(opt.paramsel.lambdas);

eps = 0;

for iter =  1:opt.Maxiter
A = newA;    
[U,S]  = eig(A);
Sigma = diag(S);
Yrotated = U'*YEncoded';
%find W ,A is fixed
for i = 1:T
	BBB = diag(D+I*(lambda/Sigma(i)));  
	BB = 1./BBB;
	B = diag(BB);
	Wthilde(i,:) = Yrotated(i,:)*X'*V*B*V';
end
W = U*Wthilde;
%find A, W is fixed
[base,value,basein] = svd(W);
ll = diag(diag(value))+eps*eye(T);
newA = base*(ll)*base';
newA = newA./trace(newA);
 end
cfr.W = (pinv(A)*W)';
cfr.A = A;
cfr.C = [];
cfr.X = [];
cfr.code = code;
end




