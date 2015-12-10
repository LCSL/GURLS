%%% set up work directory, install package
cd('/Users/Jeremiah/GitHub/GURLS/gurls/optimizers/mkl/');
% run('./utils/gurls_install.m'); savepath;
addpath('./func/'); savepath;

%%% 1. read datasets
filename = dir('./data/*.csv');
filename = {filename.name};

X = csvread(['./data/', filename{1}]);

y = X(:, 1);
X(:, 1) = [];


%%% 2. artificial example
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

%%% compute Gaussian kernels with (2\sigma^2) = 1.2^(0:49)
sigma_list = sqrt((1.2.^(0:49))/2);
K_train = rbf_array(X_train, [], sigma_list);
K_test = rbf_array(X_train, X_test, sigma_list);

K_train2 = rbf_array(X_train, [], sigma_list(1:10));
K_test2 = rbf_array(X_train, X_test, sigma_list(1:10));

K_train3 = rbf_array(X_train, [], sigma_list(20:21));
K_test3 = rbf_array(X_train, X_test, sigma_list(20:21));

%K_test = rbf_array(X_train, X_test, sigma_list);
y = y_train;

L1_cutoff = 0.0001;
L2_ratio = 0.0001;
adapt = false;

% 50 kernels
[e_list_trn, e_list_tst, A] = rls_dual_mkl_pfbs_s(...
    K_train, y_train, K_test, y_test, 0, 1e-6, 2000);

plot(1:length(e_list_trn), e_list_trn, ...
    1:length(e_list_tst), e_list_tst)

diag(A'*A)

% 2 kernels
[e_list_trn2, e_list_tst2, A] = rls_dual_mkl_pfbs_s(...
    K_train2, y_train, K_test2, y_test, 0.001, 0.001, 1000);

plot(1:length(e_list_trn2), e_list_trn2, ...
    1:length(e_list_tst2), e_list_tst2)

diag(A'*A)

% Gaussian Kernel
[e_list_trn3, e_list_tst3, A] = rls_dual_mkl_pfbs_s(...
    K_train3, y_train, K_test3, y_test, 0, 1e-6, 10000);

plot(1:length(e_list_trn3), e_list_trn3, ...
    1:length(e_list_tst3), e_list_tst3)

diag(A'*A)

%%% single kernel rbf using GURLS
options = struct('datatype','vector', 'problem', 'regression', ...
    'algorithm', 'krls', 'kernelfun', 'rbf');
options = {'datatype:vector', 'problem:regression', ...
    'algorithm:krls', 'kernelfun:rbf'};

model = gurls_train(X_train, y_train, options);
model.paramsel
ypredicted = gurls_test(model, X_test);

sum((y_test - ypredicted).^2)/sum(y_test.^2);
