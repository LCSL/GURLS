y = sign(randn(100,1));
for d=1:10
    X = randn(100,100);
    Ks(:,:,d) = X'*X;
end
[w,al,b] = simplemkl(y,Ks,100);
