function vout = paramsel_loocvdual(X,y,opt)

%	paramsel_loocvdual(X,y,opt)
%	Performs parameter selection when the dual formulation of RLS is used.
%	The leave-one-out approach is used.
%
%	NEEDS:	
%		- opt.kernel.type
%		- opt.kernel.K
%		- opt.nlambda

% Decide what you want to dump

[n,T]  = size(y);
tot = opt.nlambda;
[Q,L] = eig(opt.kernel.K);
Q = double(Q);
L = double(diag(L));
Qty = Q'*y;
filtered = L(L > 200*eps^0.5);
lmin = min(filtered)/n;
lmax = max(filtered)/n;
q = (lmax/lmin)^(1/tot);
guesses = zeros(1,tot);
LOOSQE = zeros(tot,T);
for i = 1:tot
	guesses(i) = lmin*(q^i);
	%C = rls_dual(K,y,guesses(i));
	C = rls_eigen(Q,L,Qty,guesses(i),n);
	Z = GInverseDiagonal(Q,L,guesses(i));
	opt.pred = zeros(n,T);
	for t = 1:T
		opt.pred(:,t) = y(:,t) - (C(:,t)./Z);
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
