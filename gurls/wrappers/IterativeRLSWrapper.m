function vout = IterativeRLSWrapper(X,y,opt)

% performs parameter selction and training the iterative RLS method
% specified in opt.IterRLSFun in the primal space


if ~isfield(opt,'filter');
    opt.IterRLSFun = @rls_landweberprimal;
end


% split data for parameter selection
opt.split = split_ho(X,y,opt);
tr = opt.split{1}.tr;
va = opt.split{1}.va;
Xtr = X(tr,:);
ytr = y(tr,:);
Xva = X(va,:);
yva = y(va,:);
opt.paramsel.XtX = Xtr'*Xtr;
opt.paramsel.Xty = Xtr'*ytr;

% build guesses array
if strcmp(opt.IterRLSSeriesType,'geometric');
    guesses = unique(round(opt.IterRLSMinIter.*((opt.IterRLSMaxIter/opt.IterRLSMinIter).^((1:opt.nlambda)./opt.nlambda))));
else
    guesses = opt.IterRLSMinIter:((opt.IterRLSMaxIter-opt.IterRLSMinIter)/(opt.nlambda-1)):opt.IterRLSMaxIter;
    guesses = unique(round(guesses));
end


% initialization
vout.perf = zeros(length(guesses),1);
opt_perf = 0;
prev_guess = 0;
opt.paramsel.XtXnorm = norm(opt.paramsel.XtX);

% parameter selection
for i = 1:length(guesses);
    
    opt.paramsel.lambdas = guesses(i)-prev_guess;
    opt.rls = opt.IterRLSFun(Xtr,ytr,opt);
    
    opt.pred = pred_primal(Xva,yva,opt);
    
    perf = opt.hoperf([],yva,opt);
    vout.perf(i) = mean(perf.forho);
    
    % stops when accuracy starts lowering
    new_perf = vout.perf(i);
    if new_perf>opt_perf;
        vout.opt_guess = guesses(i);
        opt_perf = new_perf;
        vout.W = opt.rls.W;
        vout.perf_opt = opt_perf;
    elseif (opt_perf-new_perf)/abs(opt_perf)>opt.IterRLSStopTol;
        break
    end

    prev_guess = guesses(i);
    if isfield(opt.rls,'WforInit')
        opt.paramsel.f0 = opt.rls.WforInit;
    else
        opt.paramsel.f0 = opt.rls.W;
    end
end

vout.guesses = guesses(1:i);
vout.perf = vout.perf(1:i);    
