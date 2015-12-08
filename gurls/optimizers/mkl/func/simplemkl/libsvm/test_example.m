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
Y= sign(randn(nPoints,1));

lambda = 1;
opt.prec = 1e-5;

tic;
[sol,b,obj] = primal_svm(0,Y,lambda,opt);
t1=toc;

tic;
[sol2,b2,obj2] = libsvm(Y,K,1/lambda,[],1);
t2=toc;

tic;
[sol3,b3,obj3] = libsvm(Y,K,1/lambda,[],0);
t3=toc;


%kern.gamma = gamma;
%kern.type = 'rbf';
%tic;
%[sol3,b3,obj3] = libsvm(Y,X,1/lambda,kern);
%t2=toc;

fprintf('%d seconds for primal, %d for libsvm\n',t1,t2);

fprintf('diff per alphas: %.04f\n',sum(abs(sol-sol3))/numel(sol));
fprintf('diff of b: %.04f\n',sum(abs(b-b3)));
fprintf('diff of objectives: %.04f\n',sum(abs(obj-obj3)));


nPoints = 500;
X = randn(nPoints,20);
X = [X;X(1:100,:)]; %repeat first 100 points
K = rbf(X,X,0.5);
Y= sign(randn(nPoints,1));
Y = [Y;Y(1:100)];
lambda =1 

[sol,b,obj] = libsvm(Y,K,1/lambda,[],1);



fprintf('test passed\n');