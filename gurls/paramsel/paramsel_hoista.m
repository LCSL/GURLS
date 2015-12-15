function vout = paramsel_hoista(X, y, opt)
% paramsel_hoinsta(X,Y,OPT)
% Performs parameter selection INSTA for elastic net is used.
% The hold-out approach is used. 
% The performance measure specified by opt.hoperf is maximized.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- split (set by the split_* routine)
%   fields with default values set through the defopt function:
%		- nlambda
%		- hoperf
%       - nholdouts
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUTS: structure with the following fields:
% -lambdas_round: cell array (opt.nholdoutsX1). For each split a cell contains the 
%       values of the regularization parameter lambda minimizing the 
%       validation error for each class.
% -perf: cell array (opt.nholdouts). For each split a cell contains a matrix 
%       with the validation error for each lambda guess and for each class
% -guesses: cell array (opt.nholdoutsX1). For each split a cell contains an 
%       array of guesses for the regularization parameter lambda
% -lambdas: mean of the optimal lambdas across splits

if isprop(opt,'paramsel')
    vout = opt.paramsel; % lets not overwrite existing parameters.
    % unless they have the same name
    
    if isfield(opt.paramsel,'perf')
        vout = rmfield(vout,'perf');
    end
    if isfield(opt.paramsel,'guesses')
        vout = rmfield(vout,'guesses');
    end
else
    opt.newprop('paramsel', struct());
end

% load in number of iterations or relative tolerance
Niter=inf;
relthre=1e-4;
if isprop(opt, 'INSTAniter')
    Niter=opt.INSTAniter;
end
if isprop(opt, 'INSTArelthre')
    relthre=opt.INSTArelthre;
end
% load in alpha for elastic net
if isprop(opt,'INSTAalpha')
    INSTAalpha=opt.INSTAalpha;
    if INSTAalpha <= 0 || INSTAalpha > 1
        error('Invalid alpha');
    end
else
    if opt.verbose
            fprintf('\t...alpha not defined. Use default value alpha=1 for LASSO\n');
            INSTAalpha=1;
    end
end


Xtytot = X'*y;

d = size(X,2);
T = size(y,2);

% Verify the parameter selection constraint is defined proper
opt.INSTASparsitylvl=ceil(opt.INSTASparsitylvl);
if (opt.INSTASAccuReq>0&&opt.INSTASparsitylvl<d)
    error('Please define parameter selection constraint correctly')
end

for nh = 1:opt.nholdouts
	if iscell(opt.split)
		tr = opt.split{nh}.tr;
		va = opt.split{nh}.va;
	else	
		tr = opt.split.tr;
		va = opt.split.va;
	end	

	n = length(tr);
	
	XtXtr = (X(tr,:))'*X(tr,:);
	Xtytr =(X(tr,:))'*y(tr,:);
	tot = opt.nlambda;


    L = 2*max(abs(Xtytot),[],2)/(n*INSTAalpha);
	guesses = n*paramsel_lambdaguesses(L, min(n,d), n, opt);
    
    if ~isprop(opt, 'rls')
        opt.newprop('rls', struct());
    end
    
	ap = zeros(tot,T);
    apt=ap;
    Sparsity = zeros(tot,1);
    for i = 1:tot
        W = rls_ista_driver( XtXtr, Xtytr, n, guesses(i),INSTAalpha,Niter,relthre,0);
        Sparsity(i) = sum(~~W);
        opt.rls.W=W;
        opt.newprop('pred', pred_primal(X(va,:),y(va,:),opt));
        opt.newprop('perf', opt.hoperf(X(va,:),y(va,:),opt));
        for t = 1:T
            ap(i,t) = opt.perf.forho(t);
        end	
    end	
    if opt.INSTASparsitylvl<d
        apt((Sparsity>opt.INSTASparsitylvl),:)=-1;
        [~, idx] = max(apt,[],1);	
    else
        Sparsity(min(ap,[],2)<opt.INSTASAccuReq)=inf;
        [~, idx] = min(Sparsity);
    end
    
    if (max(apt,[],1)<0)||isinf(min(Sparsity))
        warning('...Cannot find sych lambda that fit the requirement. Use default requirment instead')
        [~, idx] = max(ap,[],1);
    end
    vout.lambdas_round{nh} = guesses(idx);
	vout.perf{nh} = ap;
	vout.guesses{nh} = guesses;
end	
if numel(vout.lambdas_round) > 1
	lambdas = cell2mat(vout.lambdas_round');
	vout.lambdas = mean(lambdas);
else
	vout.lambdas = vout.lambdas_round{1};
end
