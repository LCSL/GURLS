function [vout] = paramsel_hoprimalr(X,y, opt)
% paramsel_hoprimalr(X,y, OPT)
% Performs parameter selection when the primal formulation of RLS is used.
% The hold-out approach is used. 
% The eigendecomposition used to compute the regularization path is computed using a randoamized method.
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
% -forho: cell array (opt.nholdoutsX1). For each split a cell contains a matrix 
%       with the validation error for each lambda guess and for each class
% -guesses: cell array (opt.nholdoutsX1). For each split a cell contains an 
%       array of guesses for the regularization parameter lambda
% -lambdas: mean of the optimal lambdas across splits

if isprop(opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
else
    opt.newprop('paramsel', struct());
end
savevars = [];

for nh = 1:opt.nholdouts
	if strcmp(class(opt.split),'cell')
		tr = opt.split{nh}.tr;
		va = opt.split{nh}.va;
	else	
		tr = opt.split.tr;
		va = opt.split.va;
	end	

	[n,d] = size(X(tr,:));
	[n,T]  = size(y(tr,:));
	
	K = X(tr,:)'*X(tr,:);
	
	k = max(round(opt.eig_percentage*d/100),1);

	[Q,L,U] = tygert_svd(K,k); % dxd matrix
	Q = double(Q);
	L = double(diag(L));
	
	QtXtY = Q'*(X(tr,:)'*y(tr,:));
	guesses = paramsel_lambdaguesses(L, k, n, opt);
	
    if ~isprop(opt, 'rls')
        opt.newprop('rls', struct());
    end
    
	tot = opt.nlambda;
	ap = zeros(tot,T);
	for i = 1:tot
		opt.rls.W = rls_eigen(Q,L,QtXtY,guesses(i),n);
		opt.newprop('pred', pred_primal(X(va,:),y(va,:),opt));
		opt.newprop('perf', opt.hoperf(X(va,:),y(va,:),opt));
		%p{i} = perf(scores,yho,{'precrec'});
		for t = 1:T
			ap(i,t) = opt.perf.forho(t);
		end	
	end	
	[dummy,idx] = max(ap,[],1);	
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
