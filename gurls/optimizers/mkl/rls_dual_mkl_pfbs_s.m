function [e_list_trn] = ...
    rls_dual_mkl_pfbs_s (K_train, y_train, L1_cutoff, L2_ratio, adapt)

% rls_dual_mkl(X, y, opt)
% Proximal Forward Backward Splitting algorithm for elastic net MKL
%
% INPUTS:
% -K: a (n x n x M) array contain M precomputed n x n kernels
% -y: response vector, n x 1 .
% -t: tunning parameter for lasso penalty
% -u: tunning parameter for ridge penalty
%
% OUTPUT: struct with the following fields:
% -W: empty matrix
% -C: matrix of coefficient vectors of dual rls estimator for each class
% -X: training samples used by the routine

M = size(K_train, 3);
n = size(K_train, 1);
iter_max = 1000;

% 1. estimate stepsize sigma = a/2 + u
eig_list = zeros(M, 1);
for m = 1:M
    eig_list(m) = eigs(K(:, :, m), 1);
end
eig_app = sum(eig_list);
u = (L2_ratio/(1 - L2_ratio)) * eig_app/2;

sigmas = eig_app/2 + u;

% 2. initialize output container, standardize data
A = zeros(n, M); % a's for M kernels
yhat = zeros(n, 1); % Ka

K = zeros(n, n, M); % standardize input kernel
K_scale = zeros(M, 1);
for m = 1:M
    K_scale(m) = svds(K_train(:, :, m), 1);
    K(:, :, m) = K_train(:, :, m)/K_scale(m);
end

y_scale = [mean(y_train) std(y_train)];  % standardize response
y = (y_train - y_scale(1))/y_scale(2);

% 3. PFBS iteration
e_list_trn = zeros(1, iter_max + 1); 
e_list_tst = zeros(1, iter_max + 1);
e_list_trn(0) = sum(y.^2);

h = waitbar(0,'PFBS in progress, PLS STAND BY..');

for iter = 2:(iter_max + 1)
    % 3.0 store previous step
    if adapt
        A_prev = A;        
        yhat_prev = yhat;
    end
    
    e_sum = 0;
    for m = 1:M
        e_sum = e_sum + ((K(:, :, m) * A(:, m) - y)/sigmas);
    end
    
    % 3.1 update a_m (n x 1)
    for m = 1:M
        % 3.2.1 GD step
        a_0 = (1 - L2_ratio) * A(:, m) - e_sum/n;
        % 3.2.2 soft-threshold
        f_norm = sqrt(a_0' * K(:,:,m) * a_0);
        A(:, m) = a_0 * max(f_norm - L1_cutoff, 0)/f_norm;
    end
    
    % 3.2. calculate prediction error
    yhat = zeros(n, 1);    
    for m = 1:M
        yhat = yhat + K(:, :, m) * A(:, m);
    end
    
    e_list_trn(iter) = sum((yhat - y).^2)/n;
    %e_list_tst(iter) = sum((yhat - y).^2)/n;
    
    % 3.3 stepsize, progress bar update, stop condition.
    if adapt
        da = A - A_prev;
        dg = (yhat - yhat_prev);
        sigmas = (n/M) * sum(da' * dg)/(dg' * dg);
    end
    
    waitbar(iter / iter_max)
        
    if (abs(e_list_trn(iter) - e_list_trn(iter-1)) < eps ||...
            e_list_trn(iter) > e_list_trn(iter-1))
        break   
    end
end

close(h)

%diag(A'*A)
%plot(1:iter_max, e_list)

end
