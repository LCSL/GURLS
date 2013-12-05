function vout = KernelIterativeRLSWrapper(X,y,opt)

% performs parameter selction and training the iterative RLS method
% specified in opt.IterRLSFun in the dual space

if ~isfield(opt,'kernelfun');
    opt.kernelfun = @kernel_linear;
end
if ~isfield(opt,'filter');
    opt.IterRLSFun = @rls_landweberdual;
end

vout = opt;

% split data for parameter selection
opt.split = split_ho(X,y,opt);
tr = opt.split{1}.tr;
va = opt.split{1}.va;
Xtr = X(tr,:);
ytr = y(tr,:);
Xva = X(va,:);
yva = y(va,:);

kernel = opt.kernelfun(X,[],opt);
opt.kernel.K = kernel.K(tr,tr);
opt.predkernel.K = kernel.K(va,tr);


% build guesses array
if strcmp(opt.IterRLSSeriesType,'geometric');
    guesses = unique(round(opt.IterRLSMinIter.*((opt.IterRLSMaxIter/opt.IterRLSMinIter).^((1:opt.nlambda)./opt.nlambda))));
else
    guesses = opt.IterRLSMinIter:((opt.IterRLSMaxIter-opt.IterRLSMinIter)/(opt.nlambda-1)):opt.IterRLSMaxIter;
    guesses = unique(round(guesses));
end

% initialization
% vout.kernelfun = opt.kernelfun;
vout.X = X(tr,:);
opt_perf = 0;
prev_guess = 0;
opt.paramsel.Knorm = norm(opt.kernel.K);

% parameter selection
for i = 1:length(guesses);
    
    opt.paramsel.lambdas = guesses(i)-prev_guess;
    opt.rls = opt.IterRLSFun(Xtr,ytr,opt);

    opt.pred = pred_dual(Xva,yva,opt);
    
    new_perf = mean(getfield(opt.hoperf([],yva,opt),'forho'));
    
    % stops when accuracy starts lowering
    if new_perf>opt_perf;
        vout.optimal_guess = guesses(i);
        opt_perf = new_perf;
        vout.C = opt.rls.C;
    elseif (opt_perf-new_perf)/abs(opt_perf)>opt.IterRLSStopTol;
        break
    end

    prev_guess = guesses(i);
    if isfield(opt.rls,'CforInit')
        opt.paramsel.f0 = opt.rls.CforInit;
    else
        opt.paramsel.f0 = opt.rls.C;
    end
end

% vout.guesses = guesses(1:i);
% vout.perf = vout.perf(1:i);    
