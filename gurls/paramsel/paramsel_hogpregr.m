function vout = paramsel_hogpregr(X,y,opt)
% paramsel_hogpregr(X,Y,OPT) 
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
%		- nlambda
%       - hoperf
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -guesses: matrix of guesses for the regularization parameter lambda
% -sigmas: array of guesses for kernel parameter sigma
% -perf: perf{1} is a matrix with the validation error for each lambda guess 
%        and for each class
% -lambdas: array of values of the regularization parameter lambda
%           minimizing the validation error for each class
% -sigma: value of the kernel parameter minimizing the validation error

if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
end

K = opt.kernel.K;

for nh = 1:opt.nholdouts

	if strcmp(class(opt.split),'cell')
		tr = opt.split{nh}.tr;
		va = opt.split{nh}.va;
	else	
		tr = opt.split.tr;
		va = opt.split.va;
	end	
	
	[n,T]  = size(y(tr,:));
    
    opt.kernel.K = K(tr,tr);
    
    opt.predkernel.K = K(va,tr);
    opt.predkernel.Ktest = diag(K(va,va));

    tot = opt.nlambda;
    
    if isfield(opt,'lambdamin')
        lmin = opt.lambdamin;
    else
        lmin = 0.001;
    end
    if isfield(opt,'lambdamax')
        lmax = opt.lambdamax;
    else
        lmax = 10;
    end
    powers = linspace(0,1,tot);
    guesses = lmin.*(lmax/lmin).^(powers);
    
    for i = 1:tot
        opt.paramsel.lambdas = guesses(i);
        opt.rls = rls_gpregr(X(tr,:),y(tr,:),opt);
        tmp = pred_gpregr(X(va,:),y(va,:),opt);
        opt.pred = tmp.means;
        opt.pred = tmp;
        
        opt.perf = opt.hoperf([],y(va,:),opt);
        for t = 1:T
            perf(i,t) = opt.perf.forho(t);
        end	
    end	
    [dummy,idx] = max(perf,[],1);	

    vout.lambdas_round{nh} = guesses(idx);
	vout.perf{nh} = perf;
	vout.guesses{nh} = guesses;
end

if numel(vout.lambdas_round) > 1
	lambdas = cell2mat(vout.lambdas_round');
	vout.lambdas = mean(lambdas);
else
	vout.lambdas = vout.lambdas_round{1};
end
