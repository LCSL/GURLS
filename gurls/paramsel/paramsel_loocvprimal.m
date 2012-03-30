function vout = paramsel_loocvprimal(X,y,opt)

%	paramsel_loocvprimal(X,y,opt)
%	Performs parameter selection when the primal formulation of RLS is used.
%	The leave-one-out approach is used.
%
%	NEEDS:	
%		- opt.nlambda

% Decide what you want to dump

K = X'*X;
[n,T]  = size(y);
d = size(X,2);
tot = opt.nlambda;
[Q,L] = eig(K);
L = double(diag(L));


tot = opt.nlambda;
guesses = paramsel_lambdaguesses(L, min(n,d), n, opt);

LEFT = X*Q;
RIGHT = Q'*X'*y;

right = Q'*X';

for i = 1:tot
	LL = L + (n*guesses(i));
	LL = LL.^(-1);
	LL = diag(LL);
	num = y - LEFT*LL*RIGHT;
	den = zeros(n,1);
	for j = 1:n
		den(j) = 1-LEFT(j,:)*LL*right(:,j);
	end	
	opt.pred = zeros(n,T);
	for t = 1:T
		opt.pred(:,t) = y(:,t) - (num(:,t)./den);
	end
	opt.perf = opt.hoperf([],y,opt);
	for t = 1:T
		ap(i,t) = opt.perf.forho(t);
	end	

end
[dummy,idx] = max(ap,[],1);	
vout.lambdas = 	guesses(idx);
vout.looe{1} = 	ap;
vout.guesses = 	guesses;
