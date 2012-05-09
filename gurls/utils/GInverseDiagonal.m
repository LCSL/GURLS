function Z = GInverseDiagonal( Q, L, lambda )

n = size(Q, 1);
t = size(lambda, 2);
Z = zeros(n, t);

D = Q.^(2);
for i = 1 : t
    d = L + (n*lambda(i));
    d  = d.^(-1);
    Z(:,i) = D*d;
end

end

