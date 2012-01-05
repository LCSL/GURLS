function vout = paramsel_loocvprimal(X,y,opt)

%	paramsel_loocvprimal(X,y,opt)
%	Performs parameter selection when the primal formulation of RLS is used.
%	The leave-one-out approach is used.
%
%	NEEDS:	
%		- opt.nlambda

% Decide what you want to dump
%savevars = {'LOOSQE','guesses'};
savevars = {};

K = X'*X;
[n,T]  = size(y);
tot = opt.nlambda;
[Q,L] = eig(K);
L = diag(L);
filtered = L(L > 200*eps^0.5);
lmin = min(filtered)/n;
lmax = max(filtered)/n;
q = (lmax/lmin)^(1/tot);

guesses = zeros(1,tot);
LOOSQE = zeros(tot,T);

LEFT = X*Q;
RIGHT = Q'*X'*y;

right = Q'*X';

for i = 1:tot
	guesses(i) = lmin*(q^i);
	LL = L + (n*guesses(i));
	LL = LL.^(-1);
	LL = diag(LL);
	num = y - LEFT*LL*RIGHT;
	%den = diag(eye(n) - LEFT*LL*right);
	den = zeros(n,1);
	for j = 1:n
		den(j) = 1-LEFT(j,:)*LL*right(:,j);
	end	
	Le = zeros(n,T);
	for t = 1:T
		Le(:,t) = num(:,t)./den;
	end
	LOOSQE(i,:) = sum(Le.*Le);
end
[dummy,bL] = min(LOOSQE);
vout.lambdas = guesses(bL);

vout.guesses = guesses;
vout.LOOSQE = LOOSQE;

% This is awesome
if numel(savevars) > 0
	[ST,I] = dbstack();
	save(ST(1).name,savevars{:});
end	
	
