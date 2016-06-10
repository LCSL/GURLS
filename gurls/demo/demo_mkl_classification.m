%%% set up work directory, install package
% run('../utils/gurls_install.m'); savepath;
addpath('./func/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MKL classification: ionosphere %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 read in ionosphere data ----
dat = csvread('./data/ionosphere.csv');
y = dat(:, 1);
X = dat(:, 2:size(dat, 2));

idx_va = randsample(size(y, 1), round(size(y, 1)/5));
idx_tr = setdiff(1:size(y, 1), idx_va);

X_tr = X(idx_tr, :);
y_tr = y(idx_tr);

X_va = X(idx_va, :);
y_va = y(idx_va);

% 2 train/test pipeline under GURLS ---- 
name = 'demo_mkl_class';
opt = gurls_defopt(name);
opt = gurls_defopt_mkl(opt);

% specify: 
% (1) hoperf => macroavg
% (2) kernel type/parameter
opt.hoperf = @perf_macroavg;
opt.mkl.type = ...
    {{'kernel_rbf', 1:0.2:4}, ...
    {'kernel_linear', 0}};

opt.seq = {...
    'split:ho', ...
    'kernel:mkl', ...    
    'paramsel:homkl', ...
    'rls:dual_mkl', ...
    'predkernel:traintest_mkl', ...
    'pred:dual_mkl', ...
    'perf:macroavg'};

opt.process{1} = [2,2,2,2,0,0,0];
opt.process{2} = [3,3,3,3,2,2,2];
gurls(X_tr, y_tr, opt, 1);
gurls(X_va, y_va, opt, 2);

% result summary/visualization

plot_mkl_path(X_tr, y_tr, opt, 'norm');
plot_mkl_path(X_tr, y_tr, opt, 'perf');

sum(opt.paramsel.norm_path(:,:,1), 1);
opt.perf