function vout = paramsel_gpregrLambdaGrid(X,y,opt)
% paramsel_gpregrLambdaGrid(X,Y,OPT) TODO to be completed
% Performs parameter selection for gaussian process regression.
% The hold-out approach is used.
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
% -looe: loo{1} is a matrix with the validation error for each lambda guess 
%        and for each class
% -guesses: array of guesses for the regularization parameter lambda 

if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
end

n = size(y,1);
        
[Q,L] = eig(opt.kernel.K);
L = double(diag(L));
	
tot = opt.nlambda;
guesses = paramsel_lambdaguesses(L, n, n, opt);
    
for i = 1:tot
    opt.paramsel.lambdas = guesses(i);
    opt.rls = rls_gpregr(X,y,opt);
    perf(i,:) = opt.rls.logprob;

end	
[dummy,idx] = max(perf,[],1);	

vout.lambdas = guesses(idx);
vout.perf = perf;
vout.guesses = guesses;

