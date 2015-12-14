function [vout] = paramsel_homkl(X,y,opt)
% paramsel_hopdual(X,y, OPT)
% Performs parameter selection on elastic net parameters for MKL
% The hold-out approach is used.
% The performance measure specified by opt.hoperf is maximized.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- split (set by the split_* routine)
%		- kernel.K_mkl (set by the kernel_mkl routine)
%   fields (mkl-related) that need to be set by hand:
%       - mkl.parrange (1x2 cell: L1/L2 penality grid)
%       or
%       - mkl.npar (1x2 int cell: num of candidate penalty terms)
%   fields (mkl-related) that need to be set by defopt_mkl:
%       - mkl.smallnumber: lower limit of non-zero L1/L2 penalty term
%       - mkl.verbose: whether to print progress bar
%       - mkl.strategy: whether to use continuation strategy in selection
%
%   fields with default values set through the defopt function:
%		- smallnumber
%		- hoperf
%       - nholdouts
%
% OUTPUTS: structure vout (opt.paramsel) with the following fields:
% -par_mkl: 2x1 cell of L1/L2 par that minimize empirical performance
% -par_path: M x tot1 x tot2 x nh double of estimated dual norm 
%            over L1/L2 grid

if isprop(opt,'paramsel')
    vout = opt.paramsel; % lets not overwrite existing parameters.
    % unless they have the same name
else
    opt.newprop('paramsel', struct());
end
vout.guesses_mkl = {};
vout.path_mkl = {};
vout.par_mkl = {};

% set the grid of L1/L2 penalty terms for elastic net
n = size(opt.kernel.K_mkl, 1);
M = size(opt.kernel.K_mkl, 3);
eig_list = opt.kernel.eig_mkl;

verbose = opt.mkl.verbose.paramsel;
iter_max = opt.mkl.iter_max.paramsel;
crit = opt.mkl.crit.paramsel;
cont_strategy = opt.mkl.strategy;

if ~isfield(opt.mkl, 'parrange')
    % estimate grid parrange if not available
    % (this code should only be executed once)
    [tot1, tot2] = opt.mkl.npar{:};
    guesses_L1 = paramsel_L1ratioguesses(y, opt, eig_list, tot1);
    guesses_L2 = [0, ...
        exp(linspace(log(opt.mkl.smallnumber), log(0.1/n), tot2))];
else
    % set grid directly using parrange
    [guesses_L1, guesses_L2] = opt.mkl.parrange{:};
    tot1 = length(guesses_L1);
    tot2 = length(guesses_L2);
end

% initiate output container for par search
ho_perf = zeros(tot1, tot2, opt.nholdouts);
ho_path = ... % norm of solution A for each kernel
    zeros(tot1, tot2, M, opt.nholdouts);

if isequal(opt.hoperf, @perf_rmsestd)
    perf_null = -1;
    perf_best = -1;
elseif isequal(opt.hoperf, @perf_macroavg)
    perf_null = 0;
    perf_best = 0;
else
    error('only ''rmsestd'' and ''macroavg'' are supported')
end

% calculate grid fit for each hold-out sample
for nh = 1:opt.nholdouts
    % specify holdout sample indices
    if iscell(opt.split)
        tr = opt.split{nh}.tr;
        va = opt.split{nh}.va;
    else
        tr = opt.split.tr;
        va = opt.split.va;
    end
    
    % prepare holdout sample/kernel
    ntr = length(tr);
    nva = length(va);
    ytr = y(tr);
    yva = y(va);
    Ktr = zeros(ntr, ntr, M);
    Kva = zeros(nva, ntr, M);
    for m = 1:M
        Ktr(:, :, m) = opt.kernel.K_mkl(tr, tr, m);
        Kva(:, :, m) = opt.kernel.K_mkl(va, tr, m);
    end
    opt.newprop('predkernel.K_mkl', Kva);
    
    %set up progress bar
    if verbose
        fprintf('\n');
        cpb = ConsoleProgressBar();
        cpb.setMinimum(0); cpb.setMaximum(tot1*tot2);
        cpb.setText(sprintf('fold%d:initializing..', nh));
        cpb.start();
        cpb_id = 0;
    end
    
    % fitting on current holdout sample over the L1/L2 penalty grid
    if ~isprop(opt, 'rls')
        opt.newprop('rls', struct());
    end
    % for each L2 term, calculate L1 solution path
    for i2 = 1:tot2
        C_mkl = zeros(ntr, M); % initiate C_mkl
        for i1 = tot1:-1:1
            % fit with/without continuation strategy
            if ~cont_strategy
                % if don't use continuation strategy
                % then reset C_mkl
                C_mkl = [];
            end
            [C_mkl] = rls_dual_mkl_pfbs(...
                Ktr, ytr, guesses_L1(i1), guesses_L2(i2), ...
                C_mkl, max(eig_list), iter_max, crit, false);
            
            % prediction
            opt.rls.C_mkl = C_mkl;
            opt.newprop('pred', pred_dual_mkl(X,yva,opt));
            opt.newprop('perf', opt.hoperf(X,yva,opt));
            
            % store result
            ho_perf(i1, i2, nh) = opt.perf.forho;
            ho_path(i1, i2, :, nh) = ...
                diag(sqrt(C_mkl'*C_mkl));
            if verbose
                %update progress bar
                cpb_id = cpb_id + 1;
                cpb.setValue(cpb_id);
                cpb.setText(...
                    sprintf('fold%d:L2(%d/%d),L1(%d/%d):%0.3f', ...
                    nh, i2, tot2, i1, tot1, opt.perf.forplot));
            end
            
            % refine L1/L2 grid
            if i2 == 1 && opt.perf.forho == perf_null
                % if null result even when result L2 = 0, break
                % fprintf('\n L1 grid reduced: %d=>%d\n', tot1, i1 - 1);
                tot1 = i1 - 1;                
            elseif i1 == 1 && opt.perf.forho == perf_null
                % if null result even when result L1 = 0, break
                % fprintf('\n L2 grid reduced: %d=>%d\n', tot2, i2 - 1);
                tot2 = i2 - 1;           
            end
        end
    end
    
    if verbose
        % destroy progress bar
        cpb.stop(); fprintf('\n');
    end
end

% refine output based on (possibly) modified tot1/tot2
guesses_L1 = guesses_L1(1:tot1);
guesses_L2 = guesses_L2(1:tot2);
ho_perf = ho_perf(1:tot1, 1:tot2, :);
ho_path = ho_path(1:tot1, 1:tot2, :, :);

% summarize perf over holdout folds then select optimal L1/L2
perf_sum = median(ho_perf, 3);
path_sum = ho_path;

[max_num, max_idx] = max(perf_sum(:));
[i1_opt, i2_opt] = ...
    ind2sub(size(perf_sum), max_idx);

% store result into paramsel.par_mkl
vout.path_mkl = reshape(...
    path_sum(1:tot1, i2_opt, :, :), tot1, M, opt.nholdouts);
vout.guesses_mkl = {guesses_L1, guesses_L2};
vout.par_mkl = {guesses_L1(i1_opt), guesses_L2(i2_opt)};
vout.cont_strategy = cont_strategy;
end
