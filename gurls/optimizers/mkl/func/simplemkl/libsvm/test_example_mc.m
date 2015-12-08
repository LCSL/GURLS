clear all

nPoints = 1500;
gamma = 1;
global K
if 1
    K = randn(nPoints,nPoints); 
    K=K'*K; 
else
    X = randn(nPoints,20);
    K = calc(kernel('rbf',gamma),data(X),data(X));
end
Y= round(10 * rand(nPoints,1));
Y = mod(Y,5) +1 ;

lambda = 1;
opt.prec = 1e-5;


tic;
[sol2,b2,obj2] = libsvm(Y,K,1/lambda,[],1);
t2=toc;

classes = unique(Y);
number_classes = numel(classes);
for i=1:number_classes
    y = Y;
    y( Y~=classes(i)) = -1;
    y( Y==classes(i)) = 1;

    [sol3(:,i),b3(i),obj3(i,1)] = libsvm(y,K,1/lambda,[],1);

    [sol(:,i),b(i),obj(i)] = primal_svm(0,y,lambda,opt);
    
end

assert(all(size(sol)==size(sol2)));
assert(all(size(sol)==size(sol3)));

assert(sum(sum(abs(sol3-sol2)))==0);
assert(sum(abs(obj3-obj2))==0);
assert(sum(abs(b3-b2))==0);

assert(mean(mean(abs(sol-sol2)))<1e-5);

fprintf('test passed\n');
