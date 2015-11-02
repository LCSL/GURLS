function [vout] = paramsel_hodual_nystrom(X,y,opt)
% paramsel_hopdual(X,y, OPT)
% Performs parameter selection when the dual formulation of RLS is used.
% The hold-out approach is used.
% The performance measure specified by opt.hoperf is maximized.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- split (set by the split_* routine)
%		- kernel.K (set by the kernel_* routines)
%   fields with default values set through the defopt function:
%		- nlambda
%		- smallnumber
%		- hoperf
%       - nholdouts
%		- kernel.type
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
vout.guesses = {};

for nh = 1:opt.nholdouts
    if iscell(opt.split)
        tr = opt.split{nh}.tr;
        va = opt.split{nh}.va;
    else
        tr = opt.split.tr;
        va = opt.split.va;
    end
    
    [n,T]  = size(y(tr,:));
    
    % TO DO: Support for linear kernel
%     if strcmp(opt.kernel.type,'linear')
%         d = size(X(tr,:),2);
%         r = min(n,d);
%     else
%         r = n;
%     end

	[ opt.kernel , nystrom ] = kernel_rbf_nystrom(X(tr,:),y(tr,:),opt);
    opt.nystrom.sampledIdx = tr(nystrom.sampledIdx);
    
    r = opt.nystrom.m;
    
    [~,L] = eig(opt.kernel.Knm' * opt.kernel.Knm);
    L = double(diag(L));
    
    if ~exist('vout','var') || ~isfield(vout, 'regrange')
        tot = opt.nlambda;
        guesses = paramsel_lambdaguesses(L, r, n, opt);
    else
        tot = numel(vout.regrange);
        guesses = vout.regrange;
    end
    
    ap = zeros(tot,T);
    Knmty = opt.kernel.Knm' * y(tr,:);
    
    if ~isprop(opt, 'rls')
        opt.newprop('rls', struct());
    end
    
    for i = 1:tot
        
        opt.rls.C = ( opt.kernel.Knm' * opt.kernel.Knm + n * guesses(i) * opt.kernel.Kmm) \ Knmty;
        
        if size(X,1) > 0
            Xva = X(va,:);
            opt.rls.X = X;
        else
            Xva = [];
        end
        
        yva = y(va,:);
        
%         if strcmp(opt.kernel.type,'linear')
%             opt.rls.W = X(tr,:)'*opt.rls.C;
%         else
%             opt.newprop('predkernel.Knm', opt.kernel.Knm(va,:));
%         end

        pk = predkernel_traintest_nystrom(Xva, yva, opt);
        opt.newprop('predkernel.Knm', pk.Knm);

        opt.newprop('pred', pred_dual_nystrom(Xva,yva,opt));
        opt.newprop('perf', opt.hoperf(Xva,yva,opt));

        for t = 1:T
            ap(i,t) = opt.perf.forho(t);
        end
    end
    [~, idx] = max(ap,[],1);
    vout.lambdas_round{nh} = guesses(idx);
    vout.perf{nh} = ap;
    vout.guesses{nh} = guesses;
end
if numel(vout.lambdas_round) > 1
    lambdas = cell2mat(vout.lambdas_round');
    vout.lambdas = median(lambdas);
else
    vout.lambdas = vout.lambdas_round{1};
end
