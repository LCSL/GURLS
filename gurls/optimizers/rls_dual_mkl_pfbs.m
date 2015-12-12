function [A, epath_trn] = ...
    rls_dual_mkl_pfbs(...
    K_train, y_train, L1_cutoff, L2_ratio, ...
    A_init, eig_app, iter_max, crit, verbose)

% rls_dual_mkl(X, y, opt)
% Proximal Forward Backward Splitting algorithm for elastic net MKL
%
% INPUTS:
% -K_train: a (n x n x M) array contain M precomputed n x n kernels
% -y_train: response vector, n x 1 .
% -L1_cutoff: tunning parameter for lasso penalty
% -L2_ratio: tunning parameter for ridge penalty
% -A_init: initialization for MKL dual solution
% -iter_max: maximum number of iteration for PFBS, default 1e4
% -crit: stopping criteria for training error, default 1e-5
% -eig_app: maximum eigenvalue for kernel matrix
%
% OUTPUT: struct with the following fields:
% -A: calculated MKL dual solution
% -epath_trn: training error over epochs

M = size(K_train, 3);
n = size(K_train, 1);

% 1.prepare data ===
K = K_train;
y = y_train;

adapt = false;

% 2. estimate stepsize sigma = a/2 + u ===

if isempty(eig_app)
    eig_list = zeros(M, 1);
    for m = 1:M
        eig_list(m) = eigs(opt.kernel.K_mkl(:, :, m), 1);
    end
    eig_app = max(eig_list);
end

sigmas = eig_app/10 + u;
u = (L2_ratio/(1 - L2_ratio)) * eig_app/2;

% 3. PFBS iteration ===
if isempty(A_init) % n x M container for parameter
    A = zeros(n, M);
else
    if any(size(A_init) ~= [n, M])
        warning(...
            'wrong A_init dimension [expect %dx%d instead of %dx%d]\n',...
            n, M, size(A_init, 1), size(A_init, 2));
        A = zeros(n, M);
    else
        A = A_init;
    end
end

if isempty(crit) % stopping criteria
    crit = 1e-5;
end

if isempty(iter_max) % stopping criteria
    iter_max = 1e4;
end

yhat = zeros(n, 1); % in-sample prediction

e_list_trn = zeros(1, iter_max + 1); % in-sample prediction error
e_list_trn(1) = sum(y.^2);

if verbose
    %set up progress bar
    fprintf('\n');
    cpb = ConsoleProgressBar();
    cpb.setMinimum(1); cpb.setMaximum(iter_max);
    cpb.setText(sprintf('L1=%0.3f,L2=%0.3f', L1_cutoff, L2_ratio));
    cpb.start();
end

for iter = 2:(iter_max + 1)
    % 3.1 update a_m (n x 1) ---
    for m = 1:M
        % 3.2.1 GD step
        a_0 = (1 - L2_ratio) * A(:, m) - (yhat - y)/(sigmas*n);
        % 3.2.2 soft thresholding
        f_norm = sqrt(a_0' * K(:,:,m) * a_0);
        A(:, m) = a_0 * max(f_norm - L1_cutoff, 0)/f_norm;
    end
    
    % 3.2. calculate train/pred error ---
    % train
    yhat = zeros(n, 1);
    for m = 1:M
        yhat = yhat + K(:, :, m) * A(:, m);
    end
    e_list_trn(iter) = sum((yhat - y).^2)/sum(y.^2);
    
    % 3.3 stepsize, progress bar update, stop condition ---
    if adapt
        da = A - A_prev;
        dg = yhat - yhat_prev;
        sigmas = (n/M) * sum(da' * dg)/(dg' * dg);
    end
    
    if verbose
        %update progress bar
        cpb.setValue(iter-1);
    end
    
    % stopping condition
    if (abs(e_list_trn(iter) - e_list_trn(iter-1)) < crit)
        break
    end
end

if verbose
    % destroy progress bar
    cpb.stop(); fprintf('\n');
end

% clean up response
if iter >= (iter_max + 1)
    warning('maximum iteration (%d) reached before convergence', iter_max)
end

epath_trn = e_list_trn(2:iter);

%diag(A'*A) % A norm for each kernel
%plot(1:length(epath_trn), epath_trn) %visualization

end
