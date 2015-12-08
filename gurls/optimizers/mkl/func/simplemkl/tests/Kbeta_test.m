
%
% Test the function Kbeta
%
clear all
K = 22;
N = 1000;

for i=1:K
    tX = randn(N,N); 
    tX = tX + tX';
    X(:,:,i) = tX;
end
w = rand(K,1);

t=cputime;
A = Kbeta(X,w,1);
t0 = cputime-t;

t=cputime;
A3 = Kbeta(X,w,0);
t2 = cputime-t;

t=cputime;
A2 = 0;
for i=1:K
    A2 = A2 + X(:,:,i) * w(i);
end
t1 = cputime-t;

fprintf('Kbeta faster fac: %.4f\n',t1/t0);
assert(sum(sum(abs(A-A2)))<1e-8);

fprintf('using symmetry faster fac: %.4f\n',t2/t0);
assert(sum(sum(abs(A3-A2)))<1e-8); 
disp('Test passed');
%
% Test the function Kbeta on non-symmetric output
%
clear all
K = 22;
N = 1000;
N2 = 2000;
for i=1:K
    tX = randn(N,N2); 
    X(:,:,i) = tX;
end
w = rand(K,1);

t=cputime;
A = Kbeta(X,w);
t0 = cputime-t;

t=cputime;
A2 = 0;
for i=1:K
    A2 = A2 + X(:,:,i) * w(i);
end
t1 = cputime-t;

fprintf('Kbeta faster fac: %.4f\n',t1/t0);
assert(sum(sum(abs(A-A2)))<1e-8);
disp('Test passed');

%
% Test the function Kbeta with beta(1,end) = 0
%
clear all
K = 22;
N = 1000;

for i=1:K
    tX = randn(N,N); 
    tX = tX + tX';
    X(:,:,i) = tX;
end
w = rand(K,1);
w(1) = 0;
w(10) = 0;
w(end) = 0;

t=cputime;
A = Kbeta(X,w);
t0 = cputime-t;

t=cputime;
A2 = 0;
for i=1:K
    A2 = A2 + X(:,:,i) * w(i);
end
t1 = cputime-t;

fprintf('Kbeta faster fac: %.4f\n',t1/t0);
assert(sum(sum(abs(A-A2)))<1e-8);
disp('Test passed');

%
% Test the function Kbeta with beta =1 and only one input dim
%
clear all
N = 1000;

tX = randn(N,N); 
tX = tX + tX';
X(:,:) = tX;

w = 1;

t=cputime;
A = Kbeta(X,w);
t0 = cputime-t;

t=cputime;
A2 = X * w;;
t1 = cputime-t;

fprintf('Kbeta faster fac: %.4f\n',t1/t0);
assert(sum(sum(abs(A-A2)))<1e-8);
disp('Test passed');
