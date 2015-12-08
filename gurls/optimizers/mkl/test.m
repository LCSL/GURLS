%%% set up work directory, install package
cd('/Users/Jeremiah/GitHub/GURLS/gurls/optimizers/mkl/');
% run('./utils/gurls_install.m'); savepath;
addpath('./func/'); savepath;

%%% 1. MNIST dataset. load data then sample 1000 obs
images_train = loadMNISTImages('./data/MNIST/train-images-idx3-ubyte');
labels_train = loadMNISTLabels('./data/MNIST/train-labels-idx1-ubyte');
images_test = loadMNISTImages('./data/MNIST/t10k-images-idx3-ubyte');
labels_test = loadMNISTLabels('./data/MNIST/t10k-labels-idx1-ubyte');

N = size(labels_train, 1);
p = 50;
n = 1000;
rng(100); idx_train = randsample(N, n);

X_train = images_train(:, idx_train)';
y_train = labels_train(idx_train);

%display_network(images(:,1:100)); % Show the first 100 images
%disp(labels(1:100));

%%% 2. artificial example
p = 50;
n = 500;
n2 = 200;

X_train = normrnd(0, 1, n, p);
X_test = normrnd(0, 1, n2, p);

beta = normrnd(0, 1, p, 1);
y_train = X_train * beta;
y_test = X_test * beta;

%%% compute Gaussian kernels with (2\sigma^2) = 1.2^(0:49)
sigma_list = sqrt((1.2.^((-10):29))/2);
K_train = rbf_array(X_train, [], sigma_list);
K_train2 = rbf_array(X_train, [], sigma_list(20:21));

%K_test = rbf_array(X_train, X_test, sigma_list);
y = y_train;

L1_cutoff = 0.01;
L2_ratio = 0.05;
adapt = false;

e_list1 = rls_dual_mkl_pfbs_s(...
    K_train, y_train, L1_cutoff, 0.1, false);
e_list2 = rls_dual_mkl_pfbs_s(...
    K_train2, y_train, L1_cutoff, 0.005, false);

M = length(e_list1);
plot(1:M, e_list1, 1:M, e_list2)


