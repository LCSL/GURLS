function vout = paramsel_loocvdual(X,y, opt)
% paramsel_loocvdual(X,y, OPT)
% Performs parameter selection when the dual formulation of RLS is used.
% The leave-one-out approach is used.
%
% INPUTS:
% -OPT: structure of options with the following fields:
%   -X: input data matrix
%   -Y: labels matrix
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

if isprop(opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
else
    opt.newprop('paramsel', struct());
end

[n,T]  = size(y);
tot = opt.nlambda;
[Q,L] = eig(opt.kernel.K);
Q = double(Q);
L = double(diag(L));
Qty = Q'*y;

if strcmp(opt.kernel.type,'linear')
    d = size(X,2);
    r = min(n,d);
else
    r = n;
end
guesses = paramsel_lambdaguesses(L, r, n, opt);

for i = 1:tot
	C = rls_eigen(Q,L,Qty,guesses(i),n);
	Z = GInverseDiagonal(Q,L,guesses(i));
	opt.newprop('pred',zeros(n,T));
	for t = 1:T
		opt.pred(:,t) = y(:,t) - (C(:,t)./Z);
	end
	opt.newprop('perf', opt.hoperf([],y,opt));
	for t = 1:T
		ap(i,t) = opt.perf.forho(t);
	end	

end	

[dummy,idx] = max(ap,[],1);	
vout.lambdas = 	guesses(idx);
vout.perf = ap;
vout.guesses = 	guesses;
