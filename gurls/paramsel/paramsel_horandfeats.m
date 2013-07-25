function [vout] = paramsel_horandfeats(X,y,opt)
% paramsel_horandfeat(X,Y,OPT)
% Performs parameter selection when the Random Features approach to RLS is used: 
%   Ali Rahimi, Ben Recht;
%   Random Features for Large-Scale Kernel Machines;
%   in Neural Information Processing Systems (NIPS) 2007.
% The hold-out approach is used. 
% The performance measure specified by opt.hoperf is maximized.
%
% INPUTS:
% -X: input data matrix
% -Y: labels matrix
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

if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
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
    
    
    if or(opt.randfeats.samplesize < 0, opt.randfeats.samplesize > n)
        ni = n;
    else 
        ni = opt.randfeats.samplesize;
    end

    [XtX,Xty,opt.rls.proj] = rp_factorize_large_real(X(tr,:)',y(tr,:)',opt.randfeats.D,'gaussian',ni);
    
	tot = opt.nlambda;
	[Q,L] = eig(XtX);
	Q = double(Q);
	L = double(diag(L));
	QtXtY = Q'*Xty;
    
	guesses = paramsel_lambdaguesses(L, min(n,opt.randfeats.D*2), n, opt);

	perf = zeros(tot,T);
	for i = 1:tot
		opt.rls.W = rls_eigen(Q,L,QtXtY,guesses(i),n);
		opt.pred = pred_randfeats(X(va,:),y(va,:),opt);
		opt.perf = opt.hoperf(X(va,:),y(va,:),opt);
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
