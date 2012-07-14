function vout = paramsel_loocvprimal(X,y,opt)
% paramsel_loocvprimal(X,Y,OPT)
% Performs parameter selection when the primal formulation of RLS is used.
% The leave-one-out approach is used.
%
% INPUTS:
% -X: input data matrix
% -Y: labels matrix
% -OPT: structure of options with the following fields with default values
% set through the defopt function:
%		- nlambda
%		- smallnumber
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUTS: structURE with the following fields:
% -lambdas: array of values of the regularization parameter lambda
%           minimizing the validation error for each class
% -perf: is a matrix with the validation error for each lambda guess 
%        and for each class
% -guesses: array of guesses for the regularization parameter lambda 

if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
end

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
vout.perf = 	ap;
vout.guesses = 	guesses;
