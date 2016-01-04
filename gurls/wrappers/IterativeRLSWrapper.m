function vout = IterativeRLSWrapper(X, y, opt)
% 
% IterativeRLSWrapper(X, y, opt) 
% Performs parameter selection and training the iterative RLS method
% specified in opt.IterRLS.fun in the primal space.
%
% INPUTS:
%   -X: input data matrix
%   -OPT: structure of options with the following fields (and subfields):
%       -iterRLS: options for Iterative RLS     
%
% OUTPUT:
% -vout: a struct with the optimal num of iterations and vector of parameters
%
% NOTE: This is function is probably obsolete and should be removed.

%% Initialize fields (that might be missing from opt)
if ~isa(opt, 'GurlsOptions')
    opt = GurlsOptions(opt);
end
if ~isprop(opt, 'split')
    opt.newprop('split', struct());
end
if ~isprop(opt, 'paramsel')
    opt.newprop('paramsel', struct());
end
if ~isprop(opt, 'rls')
    opt.newprop('rls', struct());
end
if ~isprop(opt, 'pred')
    opt.newprop('pred', struct());
end

if ~isprop(opt, 'filter')
    opt.iterRLS.fun = @rls_landweberprimal;
end

%% split data for parameter selection
opt.split = split_ho(X,y,opt);
tr = opt.split{1}.tr;
va = opt.split{1}.va;
Xtr = X(tr,:);
ytr = y(tr,:);
Xva = X(va,:);
yva = y(va,:);
opt.paramsel.XtX = Xtr'*Xtr;
opt.paramsel.Xty = Xtr'*ytr;

%% output struct
% make a new empty struct
vout = struct; 
% vout = gurls_defopt(opt.name);
% vout.newprops(opt); 
% vout.iterRLS = opt.iterRLS; 

%% build guesses array
if strcmp(opt.iterRLS.seriestype, 'geometric');
    guesses = unique(round(opt.iterRLS.miniter.*((opt.iterRLS.maxiter/opt.iterRLS.miniter).^((1:opt.nlambda)./opt.nlambda))));
else
    guesses = opt.iterRLS.miniter:((opt.iterRLS.maxiter - opt.iterRLS.miniter)/(opt.nlambda-1)):opt.iterRLS.maxiter;
    guesses = unique(round(guesses));
end

%% initialization
opt_perf = 0;
% prev_guess = 0;
opt.paramsel.XtXnorm = norm(opt.paramsel.XtX);

%% parameter selection
for iter = 1:length(guesses);
    
    opt.paramsel.lambdas = guesses(iter); % - prev_guess;
    opt.rls = opt.iterRLS.fun(Xtr, ytr, opt);
    opt.pred = pred_primal(Xva, yva, opt);
    
    new_perf = mean(getfield(opt.hoperf([], yva, opt), 'forho'));
    % perf = opt.hoperf([], yva, opt);
    % new_perf = mean(perf.forho);
    
    %% stops when accuracy starts lowering
    % fprintf('%d, %d, %f\n', iter, guesses(iter), new_perf);
    if new_perf>opt_perf
        
        opt_perf = new_perf;
        vout.optimal_guess = guesses(iter);
        vout.W = opt.rls.W;
        % vout.perf_opt = opt_perf;
        
    elseif (opt_perf - new_perf)/abs(opt_perf) > opt.iterRLS.stoptol
        break
    end
    
    % prev_guess = guesses(iter);
    if isfield(opt.rls, 'WforInit')
        opt.paramsel.f0 = opt.rls.WforInit;
    else
        opt.paramsel.f0 = opt.rls.W;
    end
    
end

% vout.guesses = guesses(1:i);
% vout.perf = vout.perf(1:i);
