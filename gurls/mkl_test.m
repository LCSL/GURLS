%%% set up work directory, install package
cd('/Users/Jeremiah/GitHub/GURLS/gurls/');
% run('./utils/gurls_install.m'); savepath;
addpath('./func/'); savepath;

%%%% 0. data generation ====
p = 50;
n = 500;
n2 = 200;

X = normrnd(0, 1, n + n2, p);
beta = normrnd(0, 1, p, 1);
y = X * beta;

X_train = X(1:n, :);
y_train = y(1:n);

X_test = X((n+1):(n+n2), :);
y_test = y((n+1):(n+n2));

%%%% 1. process ====

% initialization
opt = gurls_defopt('mkl_test');
opt.newprop('mkl', struct());

% kernel step (add opt.type)
opt.mkl.type = ...
    {{'kernel_rbf', sqrt((1.2.^(0:49))/2)}, {'kernel_linear', 0}};
opt.kernel = kernel_mkl(X, y, opt); 

% split step (change opt.nholdouts 1 -> 10)
opt.nholdouts = 10;
opt.newprop('split', split_ho(X, y, opt)); 

% paramsel step (add opt.mkl.L1range and opt.mkl.L2range)
opt.mkl.L1range
opt.mkl.L2range 



