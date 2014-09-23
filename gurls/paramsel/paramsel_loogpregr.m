function vout = paramsel_loogpregr(X,y, opt)
% paramsel_loogpregr(X,y, OPT)
% Performs parameter selection for gaussian process regression.
% The leave-one-out approach is used.
%
% INPUTS:
% -OPT: structure of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- kernel.K (set by the kernel_* routines)
%   fields with default values set through the defopt function:
%		- nlambda
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -guesses: array of guesses for the regularization parameter lambda 
% -perf: perf is a matrix with the validation error for each lambda guess 
%        and for each class
% -lambdas: array of values of the regularization parameter lambda
%           minimizing the validation error for each class

if isprop(opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
else
    opt.newprop('paramsel', struct());
end    
    
[n,T]  = size(y);
tot = opt.nlambda;
K = opt.kernel.K;


if isprop(opt,'lambdamin')
    lmin = opt.lambdamin;
else
    lmin = 0.001;
end
if isprop(opt,'lambdamax')
    lmax = opt.lambdamax;
else
    lmax = 10;
end
powers = linspace(0,1,tot);
guesses = lmin.*(lmax/lmin).^(powers);


perf = zeros(tot,T);
for k = 1:n;
    tr = setdiff(1:n,k);
    opt.kernel.K = K(tr,tr);
    if ~isprop(opt, 'predkernel')
        opt.newprop('predkernel', struct());
    end
    opt.predkernel.K = K(k,tr);
    opt.predkernel.Ktest = K(k,k);
    for i = 1:tot
        opt.paramsel.lambdas = guesses(i);
        opt.newprop('rls', rls_gpregr(X(tr,:),y(tr,:),opt));
        tmp = pred_gpregr(X(k,:),y(k,:),opt);
        opt.newprop('pred',tmp.means);
        opt.newprop('perf', opt.hoperf([],y(k,:),opt));
        for t = 1:T
            perf(i,t) = opt.perf.forho(t)./n+perf(i,t);
        end	
    end
end	

opt.kernel.K = K;

[dummy,idx] = max(perf,[],1);	
vout.lambdas = 	guesses(idx);
vout.perf = 	perf;
vout.guesses = 	guesses;
