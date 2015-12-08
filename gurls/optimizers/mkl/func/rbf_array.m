function K = rbf_array(X, X2, sigma_list)
% calculate normalized rbf kernels based on an array of sigma
n = size(X, 1);
M = length(sigma_list);
K = zeros(n, n, M);

for m = 1:M
    K(:, :, m) = rbf2(X, X2, sigma_list(m));
end

end