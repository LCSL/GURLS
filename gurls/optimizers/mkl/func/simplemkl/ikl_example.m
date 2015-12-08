%
% An example script to introduce into the infinite kernel learning algorithm
%

clear all

fprintf('generating training and test data...');
nTiles = [2;4];
nPoints = 300;
noisedims = 18;

% generate data
[X,Y] = checkerboard(nTiles,nPoints,[0.5;0.5]);
[Xtest,Ytest] = checkerboard(nTiles,nPoints);

X(:,size(X,2)+(1:noisedims)) = randn(size(X,1),noisedims);
Xtest(:,size(Xtest,2)+(1:noisedims)) = randn(size(Xtest,1),noisedims);

% normalize data
X = X - ones(size(X,1),1)*mean(X);
X = X * diag(1./std(X));

Xtest = Xtest - ones(size(Xtest,1),1)*mean(Xtest);
Xtest = Xtest * diag(1./std(Xtest));

% plot the data
figure(3); clf
subplot(221)
plot(X(Y==1,1),X(Y==1,2),'ro','MarkerSize',5);
hold on
plot(X(Y==-1,1),X(Y==-1,2),'b*','MarkerSize',5);
title('chessboard training data','FontSize',20);
axis([-2,2,-2,2]); 
axis square
xlabel('Dimension 1','FontSize',16);
ylabel('Dimension 2','FontSize',16);
set(gca,'FontSize',16);
subplot(222)
plot(X(Y==1,1),X(Y==1,3),'ro','MarkerSize',5);
hold on
plot(X(Y==-1,1),X(Y==-1,3),'b*','MarkerSize',5);
axis([-2,2,-2,2]); 
axis square
xlabel('Dimension 1','FontSize',16);
ylabel('Dimension 3 (noise)','FontSize',16);
set(gca,'FontSize',16);
subplot(223)
plot(X(Y==1,3),X(Y==1,2),'ro','MarkerSize',5);
hold on
plot(X(Y==-1,3),X(Y==-1,2),'b*','MarkerSize',5);
axis([-2,2,-2,2]); 
axis square
ylabel('Dimension 2','FontSize',16);
xlabel('Dimension 3 (noise)','FontSize',16);
set(gca,'FontSize',16);
subplot(224)
plot(X(Y==1,3),X(Y==1,4),'ro','MarkerSize',5);
hold on
plot(X(Y==-1,3),X(Y==-1,4),'b*','MarkerSize',5);
axis([-2,2,-2,2]); 
axis square
xlabel('Dimension 3 (noise)','FontSize',16);
ylabel('Dimension 4 (noise)','FontSize',16);
set(gca,'FontSize',16);


%
% Start of IKL algorithm
%

% first compute the distance in each dimension separately
D = zeros(size(X,1),size(X,1),size(X,2));
for k=1:size(X,2)
    D(:,:,k) = dist_euclid(X(:,k),X(:,k));
    Dtest(:,:,k) = dist_euclid(Xtest(:,k),X(:,k));
end

% and wait for user to acknowledge that
fprintf('press key\n');
pause


fac_range = [0;10];

opts = [];
subprob_opts = [];
subprob_opts.constraints_per_iteration = 1;
subprob_opts.do_plots = 1;
subprob_opts.do_continuation = 0;
subprob_opts.do_continuation_inverse = 0;
subprob_opts.scale_factors = 1;


C = 10; 

% define the subproblem function ...
% the variables D, fac_range and subprob_opts are only defined
% here. The IKL function IKLmain will never touch them
subprob = @(method,al,lam,parm,w) IKLsubproblem_nonnegisotropic(method,al,lam,D,parm,w,fac_range,subprob_opts);
%subprob = @(method,al,lam,parm,w) IKLsubproblem_singledimgaussian(method,al,lam,D,parm,w,fac_range,subprob_opts);

% ... now train the model
[w,alpha,b,Theta] = IKLmain(Y,subprob,C,opts);

% build the kernel for test data ...
Ktest = IKLsubproblem_nonnegisotropic('build',[],[],Dtest,Theta,w);
% ... and test
Ytest_pred = sign(Ktest * alpha + b);

% print some statistics
fprintf('Algorithm selected %d kernels\n',numel(w));
fprintf('test accuracy: %.5g\n',mean(Ytest==Ytest_pred));
[ignore,ind]= max(w);
fprintf('parameters of most active kernel (w=%.2g) = ',w(ind));
fprintf(' %.2f',Theta(ind,:));
fprintf('\n');


fprintf('press key to exit\n');
pause