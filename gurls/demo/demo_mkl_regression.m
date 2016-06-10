%%% set up work directory, install package
% run('../utils/gurls_install.m'); savepath;
addpath('./func/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MKL regression on linear data %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 generate linear data----
p = 50;
n = 500;
n2 = 200;

X_0 = normrnd(0, 1, n + n2, p);
beta = normrnd(0, 1, p, 1);

sigma = 20;
K_0 = exp(-square_distance(X_0', X_0')/(sigma^2));
alpha = normrnd(0, 1, n + n2, 1);

y_0 = X_0 * beta + K_0 * alpha;

X_tr = X_0(1:n, :);
y_tr = y_0(1:n);

X_va = X_0((n+1):(n+n2), :);
y_va = y_0((n+1):(n+n2));

% 2 train/test pipeline under GURLS ----
name = 'demo_mkl_reg';
opt = gurls_defopt(name);
opt = gurls_defopt_mkl(opt);
% specify kernel type/parameter
% (sigma for rbf, none for linear) 
opt.mkl.type = ...
    {{'kernel_rbf', sqrt((1.2.^(0:49))/2)}};
opt.mkl.strategy = false;

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
gurls(X_va, y_va, opt, 2);

% result summary/visualization  
plot_mkl_path(X_tr, y_tr, opt, 'norm');
plot_mkl_path(X_tr, y_tr, opt, 'perf');
% kernel
norm_summary = median(opt.paramsel.norm_path, 3);
kernel_importance = sum(norm_summary, 1);
[norm_value, norm_order] = ...
    sort(kernel_importance, 'descend');
[norm_order(1:5); norm_value(1:5)];

% 3 train/test using true model  ----
name = 'demo_mkl_reg_true';
opt = gurls_defopt(name);
opt = gurls_defopt_mkl(opt);
% specify kernel type/parameter
% (sigma for rbf, none for linear) 
opt.mkl.type = ...
    {{'kernel_rbf', 10}, {'kernel_linear', 0}};
opt.mkl.par_mkl = {[0], [0]};

% skip paramsel
opt.seq = {...
    'split:ho', ...
    'kernel:mkl', ...    
    'rls:dual_mkl', ...
    'predkernel:traintest_mkl', ...
    'pred:dual_mkl', ...
    'perf:rmsestd'};
opt.process{1} = [2,2,2,0,0,0];
opt.process{2} = [3,3,3,2,2,2];
gurls(X_tr, y_tr, opt, 1);
gurls(X_va, y_va, opt, 2);