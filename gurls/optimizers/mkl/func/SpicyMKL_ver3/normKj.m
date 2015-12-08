    function [v,wj] = normKj(K,u,activeset)

M = size(u,2);
N = size(u,1);

if nargin < 3
    b_active = 0;
else
    M = length(activeset);
    b_active = 1;
end;


v = zeros(N,M);
wj = zeros(M,1);
if b_active
    mcount = 0;
    for m = activeset
        mcount = mcount + 1;
        v(:,mcount) = K(:,:,m)*u(:,m);
        wj(mcount) = sqrt(max(v(:,mcount)'*u(:,m),0));
    end;
else
    for m = 1:M
        v(:,m) = K(:,:,m)*u(:,m);
        wj(m) = sqrt(max(v(:,m)'*u(:,m),0));
    end;
end;