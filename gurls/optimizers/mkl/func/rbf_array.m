function K = rbf_array(X, X2, sigma_list)
% calculate normalized rbf kernels based on an array of sigma
n = size(X, 1);
if isempty(X2)
    n2 = n;
else
    n2 = size(X2, 1);
end

M = length(sigma_list);
K = zeros(n2, n, M);

for m = 1:M
    K(:, :, m) = rbf2(X, X2, sigma_list(m));
end

end