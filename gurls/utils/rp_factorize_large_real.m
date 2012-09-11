% Sample from the fourier bases of the a kernel to obtain a rank d
% approximation for the kernel.
% G is a d by N complex matrix so that G'G approximates K.
% Returns G*G' and G*y rather than just G in order to economize
% storage.
function [GG,Gy,W] = rp_factorize_large_real(X,y,d,kernel,psize)
D = size(X,1);
N = size(X,2);
T = size(y,1);

W = rp_projections(D,d,kernel);

% G = exp(i*W*X);
% But G is too big to store in memory.
% Instead, since G= [G1 ... Gn] = [exp(i*W*X1) ... exp(i*W*Xn)]
% GG' = [G1G1' + ... + GnGn'].
% Gy = G1*y1 + ... Gn*yn.
GG = zeros(d*2,d*2);
Gy = zeros(d*2,T);
for i=1:psize:N
    bend = min(i+psize-1,N);
    Gi = rp_apply_real(X(:,i:bend),W);
    yi = y(:,i:bend);
    GG = GG + Gi*Gi';
    Gy = Gy + Gi*yi';
end
end
