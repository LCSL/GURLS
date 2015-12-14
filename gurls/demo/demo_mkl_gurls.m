%%% set up work directory, install package
cd('/Users/Jeremiah/GitHub/GURLS/gurls/demo');
% run('../utils/gurls_install.m'); savepath;

%%%% 1. regression on linear data ====
% 1.1 generate linear data----
p = 50;
n = 500;
n2 = 200;

X_0 = normrnd(0, 1, n + n2, p);
beta = normrnd(0, 1, p, 1);
y_0 = X_0 * beta;

X_tr = X_0(1:n, :);
y_tr = y_0(1:n);

X_va = X_0((n+1):(n+n2), :);
y_va = y_0((n+1):(n+n2));

% 1.2 naive train/test under GURLS ----
name = 'demo_mkl_reg';
opt = gurls_defopt(name);
opt = gurls_defopt_mkl(opt);
% specify kernel type/parameter
% (sigma for rbf, none for linear) 
opt.mkl.type = ...
    {{'kernel_rbf', sqrt((1.2.^(0:49))/2)}, ...
    {'kernel_linear', 0}};

opt.seq = {...
    'split:ho', ...
    'kernel:mkl', ...    
    'paramsel:homkl', ...
    'rls:dual_mkl', ...
    'predkernel:traintest_mkl', ...
    'pred:dual_mkl', ...
    'perf:rmsestd'};

opt.process{1} = [2,2,2,2,0,0,0];
opt.process{2} = [3,3,3,3,2,2,2];
gurls(X_tr, y_tr, opt, 1);
gurls (X_va, y_va, opt, 2);

% 1.3 additional justification

% determine important kernel using 
% L1 regularization path
plot_mkl_L1path(X_tr, y_tr, opt);
sum(opt.paramsel.path_mkl(:,:,1), 1);
gurls.perf

% reset L1 parameter then re-select model
opt.mkl.parrange = {[1e-4 * (4:0.5:6)], [0]};
opt.mkl.strategy = false; %(don't use cont strategy)
save(opt.savefile, 'opt', '-v7.3');

opt.process{3} = [3,3,2,2,0,0,0];
gurls (X_tr, y_tr, opt, 3);
gurls (X_va, y_va, opt, 2);
opt.perf % slightly better performance 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2. classification: ionosphere %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1 read in ionosphere data ----
dat = csvread('./data/ionosphere.csv');
y = dat(:, 1);
X = dat(:, 2:size(dat, 2));

idx_va = randsample(size(y, 1), round(size(y, 1)/5));
idx_tr = setdiff(1:size(y, 1), idx_va);

X_tr = X(idx_tr, :);
y_tr = y(idx_tr);

X_va = X(idx_va, :);
y_va = y(idx_va);

% 2.2 naive train/test under GURLS ---- 
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

% 2.3 additional adjustment ----
% determine important kernel using 
% L1 regularization path
plot_mkl_L1path(X_tr, y_tr, opt);
sum(opt.paramsel.path_mkl(:,:,1), 1);
opt.perf

% reset L1 parameter then re-select model
opt.mkl.parrange = {[1e-4 * (8:15)], [0]};
opt.mkl.strategy = false; %(don't use cont strategy)

save(opt.savefile, 'opt', '-v7.3');

opt.process{3} = [3,3,2,2,0,0,0];
gurls (X_tr, y_tr, opt, 3);
gurls (X_va, y_va, opt, 2);
opt.perf % better performance 

