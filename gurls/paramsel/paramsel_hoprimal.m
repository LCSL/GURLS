function [vout] = paramsel_hoprimal(X, y, opt)
% paramsel_hoprimal(X,Y,OPT)
% Performs parameter selection when the primal formulation of RLS is used.
% The hold-out approach is used. 
% The performance measure specified by opt.hoperf is maximized.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- split (set by the split_* routine)
%   fields with default values set through the defopt function:
%		- nlambda
%		- smallnumber
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


%verify if matrix XtX has already been computed (especially for online RLS)
if isprop(opt,'kernel');
    if isfield(opt.kernel,'XtX');
        Ktot = opt.kernel.XtX;
    else
        Ktot = X'*X;
    end
    if isfield(opt.kernel,'XtX');
        Xtytot = opt.kernel.Xty;
    else
        Xtytot = X'*y;
    end
else
    Ktot = X'*X;
    Xtytot = X'*y;
    opt.newprop('kernel',struct());
end

d = size(X,2);
T = size(y,2);

for nh = 1:opt.nholdouts
	if iscell(opt.split)
		tr = opt.split{nh}.tr;
		va = opt.split{nh}.va;
	else	
		tr = opt.split.tr;
		va = opt.split.va;
	end	

	n = length(tr);
	
	K = Ktot - X(va,:)'*X(va,:);
	
	tot = opt.nlambda;
	[Q,L] = eig(K);
	Q = double(Q);
	L = double(diag(L));
	QtXtY = Q'*(Xtytot - X(va,:)'*y(va,:));
	
	guesses = paramsel_lambdaguesses(L, min(n,d), n, opt);
    
    if ~isprop(opt, 'rls')
        opt.newprop('rls', struct());
    end
    
	ap = zeros(tot,T);
	for i = 1:tot
		opt.rls.W = rls_eigen(Q,L,QtXtY,guesses(i),n);
		opt.newprop('pred', pred_primal(X(va,:),y(va,:),opt));
		opt.newprop('perf', opt.hoperf(X(va,:),y(va,:),opt));
		for t = 1:T
			ap(i,t) = opt.perf.forho(t);
		end	
	end	
	[dummy, idx] = max(ap,[],1);	
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
