function vout = paramsel_hogpregr(X,y,opt)
% paramsel_hogpregr(X,Y,OPT) TODO to be completed
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

tot = opt.nlambda;
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
    
%     lmax = mean(std(y(tr)));
%     lmin = mean(std(y(tr)))*10^-5;
%     guesses = lmin.*(lmax/lmin).^linspace(0,1,tot);
    
    opt.predkernel.K = K(va,tr);
    opt.predkernel.Ktest = diag(K(va,va));

    [Q,L] = eig(opt.kernel.K(tr,tr));
	L = double(diag(L));
	
	tot = opt.nlambda;
	guesses = paramsel_lambdaguesses(L, n, n, opt);
    
    
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
