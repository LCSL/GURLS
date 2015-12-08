% K = rbf(X1,X2,sig)
%
% computes:
% K = exp(-||X1-X2||/(2*sig^2));
function K = rbf(X1,X2,sig)

if numel(sig)~=1
    error('third argument (sigma) must be a scalar');
end

K = exp(-dist_euclid(X1,X2)/(2*sig^2));
