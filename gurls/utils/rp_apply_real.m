function G = rp_apply_real(X,W)
V = W*X;
G = [cos(V);
    sin(V)];
% G = X;
end