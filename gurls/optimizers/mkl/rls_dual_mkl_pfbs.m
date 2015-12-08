function [cfr] = rls_dual_mkl_pfbs (K, y, u, t)

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

M = size(K, 3);
n = size(K, 1);
iter_max = 1000;

% 1. choose stepsize (a + b)/2 + u ???????????????
% currently: different step size for each group
sigmas = zeros(1, M);
for i = 1:M
    sigmas(i) = eigs(K(:, :, i), 1)/2 + u;
end


% 2. initialize output container
A = zeros(n, M); % a's for M kernels

% 3. PFBS iteration
e_list = zeros(1, iter_max);
h = waitbar(0,'PFBS in progress, PLS STAND BY..');
for iter = 1:iter_max
    % 3.1 update error e = Ka - y (n x 1)
    e_sum = zeros(n, 1);
    
    for m = 1:M
        e_sum = e_sum + ((K(:, :, m) * A(:, m) - y)/sigmas(m));
    end
    
    % 3.2 update a_m (n x 1)
    for m = 1:M
        % 3.2.1 GD step
        a_0 = (1 - u/sigmas(m)) * A(:, m) - e_sum/n;
        % 3.2.2 soft-threshold
        f_norm = sqrt(a_0' * K(:,:,m) * a_0);
        A(:, m) = a_0 * max(f_norm - t/sigmas(m), 0)/f_norm;
    end
    
    % 3.3. calculate prediction error
    yhat = zeros(n, 1);    
    for m = 1:M
        yhat = yhat + K(:, :, m) * A(:, m);
    end
    
    e_list(iter) = sum((yhat - y).^2);
    
    % 3.4 clean up, progress bar update
    if adapt
        %u = L1_cutoff * sigmas;
        da = A - A_prev;
        dg = (yhat - yhat_prev);
        %inner_aa = sum(sum(da.^2)) * u;
        %inner_ag = sum(da' * dg) / n;
        %inner_gg = (dg' * dg * M) / (n^2);
        %sigmas = (inner_aa + inner_ag)/...
        %    (u*inner_aa + inner_gg + 2 * u * inner_ag) ;
        sigmas = (n/M) * sum(da' * dg)/(dg' * dg);
    end
    waitbar(iter / iter_max)
end
close(h)

end
