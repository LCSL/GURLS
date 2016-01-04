function vout = KernelIterativeRLSWrapper(X, y, opt)
%
% KernelIterativeRLSWrapper(X, y, opt)
% Performs parameter selection and training the iterative RLS method
% specified in opt.IterRLS.fun in the dual space.
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

if ~isprop(opt, 'kernelfun')
    if ~isfield(opt.kernel, 'type')
        opt.kernel.type = 'linear';
    end
    kernelfun = str2func(['kernel_' opt.kernel.type]);
else
    % compatibility
    % opt.kernelfun = @kernel_linear;
    kernelfun = opt.kernelfun;
end

if ~isprop(opt, 'filter')
    opt.iterRLS.fun = @rls_landweberdual;
end

%% split data for parameter selection
opt.split = split_ho(X,y,opt);
tr = opt.split{1}.tr;
va = opt.split{1}.va;
Xtr = X(tr,:);
ytr = y(tr,:);
Xva = X(va,:);
yva = y(va,:);

%%  if kernel is rbf and sigma is not provided do CV for single lambda
if strcmp(opt.kernel.type, 'rbf') && ~isfield(opt.paramsel, 'sigma');
    opt.nlambda = 1;
    opt.nsigma = 10;
    paramsel = paramsel_siglamho(X, y, opt);
    opt.paramsel.sigma = paramsel.sigma;
end

kernel = kernelfun(X, [], opt);
opt.kernel = kernel;
opt.kernel.K = kernel.K(tr, tr);
opt.predkernel.K = kernel.K(va, tr);

%% output struct
vout = struct; % opt; % gurls_defopt(opt.name);
% vout.newprops(struct('optimal_guess',[], 'C', [], 'kernelfun', kernelfun))

% vout.kernelfun = kernelfun; % keep track of kernel function (outdated),
% keep track of X, kernel (and sigma) in output
vout.X = X(tr,:);
vout.kernel = opt.kernel;
if isfield(opt.paramsel, 'sigma')
    vout.paramsel.sigma = opt.paramsel.sigma;
end

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
opt.paramsel.Knorm = norm(opt.kernel.K);

%% parameter selection
for iter = 1:length(guesses);
    
    opt.paramsel.lambdas = guesses(iter); % - prev_guess;
    opt.rls = opt.iterRLS.fun(Xtr, ytr, opt);
    opt.pred = pred_dual(Xva, yva, opt);
    
    new_perf = mean(getfield(opt.hoperf([], yva, opt), 'forho'));
    
    %% stops when accuracy starts lowering
    % fprintf('%d, %d, %f\n', iter, guesses(iter), new_perf);
    if new_perf>opt_perf
        
        opt_perf = new_perf;
        
        vout.optimal_guess = guesses(iter);
        vout.C = opt.rls.C;
        
    elseif (opt_perf - new_perf)/abs(opt_perf) > opt.iterRLS.stoptol
        break
    end
    
    % prev_guess = guesses(iter);
    if isfield(opt.rls, 'CforInit')
        opt.paramsel.f0 = opt.rls.CforInit;
    else
        opt.paramsel.f0 = opt.rls.C;
    end
    
end

% vout.guesses = guesses(1:i);
% vout.perf = vout.perf(1:i);
