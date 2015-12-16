%%% set up work directory, install package
% run('../utils/gurls_install.m'); savepath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MKL regression on linear data %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 generate linear data----
p = 50;
n = 500;
n2 = 200;

X_0 = normrnd(0, 1, n + n2, p);
beta = normrnd(0, 1, p, 1);

sigma = 10;
K_0 = exp(-square_distance(X_0', X_0')/(2*sigma^2));
alpha = normrnd(0, 1, n + n2, 1);

y_0 = X_0 * beta + K_0 * alpha;

X_tr = X_0(1:n, :);
y_tr = y_0(1:n);

X_va = X_0((n+1):(n+n2), :);
y_va = y_0((n+1):(n+n2));

% 2 naive train/test under GURLS ----
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
gurls(X_va, y_va, opt, 2);

plot_mkl_path(X_tr, y_tr, opt, 'norm');
plot_mkl_path(X_tr, y_tr, opt, 'perf');
opt.perf

% 2 naive train/test under GURLS ----
name = 'demo_mkl_reg_noStrategy';
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
gurls(X_va, y_va, opt, 2);

% 3 additional justification ----

opt.mkl.strategy = false; %(don't use cont strategy)
opt.process{1} = [2,2,2,2,0,0,0];
opt.process{2} = [3,3,3,3,2,2,2];
gurls(X_tr, y_tr, opt, 1);
gurls(X_va, y_va, opt, 2);

% reset L1 parameter then re-select model
opt.mkl.parrange = {[1e-4 * (4:0.5:6)], [0]};
opt.mkl.strategy = false; %(don't use cont strategy)
save(opt.savefile, 'opt', '-v7.3');

opt.process{3} = [3,3,2,2,0,0,0];
gurls (X_tr, y_tr, opt, 3);
gurls (X_va, y_va, opt, 2);
opt.perf % slightly better performance 

% 4 re-fit using optimal kernel ----
name = 'demo_mkl_reg_linear';
opt = gurls_defopt(name);
opt = gurls_defopt_mkl(opt);
% specify kernel type/parameter
% (sigma for rbf, none for linear) 
opt.mkl.type = {{'kernel_linear', 0}};
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