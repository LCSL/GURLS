function vout = paramsel_loocvdual(X,y,opt)
% paramsel_loocvdual(X,Y,OPT)
% Performs parameter selection when the dual formulation of RLS is used.
% The leave-one-out approach is used.
%
% INPUTS:
% -X: input data matrix
% -Y: labels matrix
% -OPT: structure of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- kernel.K (set by the kernel_* routines)
%   fields with default values set through the defopt function:
%		- kernel.type
%		- nlambda
%       - hoperf
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -lambdas: array of values of the regularization parameter lambda
%           minimizing the validation error for each class
% -perf: is a matrix with the validation error for each lambda guess 
%        and for each class
% -guesses: array of guesses for the regularization parameter lambda 

if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
end

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
vout.perf = 	ap;
vout.guesses = 	guesses;
