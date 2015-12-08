%
% An example script to introduce into the infinite kernel learning algorithm
%

clear all

% first generate some data ...
fprintf('generating training and test data...');
nTiles = [2;2];
nPoints = 300;
nTestPoints = 500;
noisedims = 5;

% ... checkerboard data ...
[X,Y] = checkerboard(nTiles,nPoints,[1;1]);
[Xtest,Ytest] = checkerboard(nTiles,nTestPoints);

% ... add noise ...
if noisedims >0
    X(:,size(X,2)+(1:noisedims)) = randn(size(X,1),noisedims);
    Xtest(:,size(Xtest,2)+(1:noisedims)) = randn(size(Xtest,1),noisedims);
end


% ... normalize data ...
X = X - ones(size(X,1),1)*mean(X);
X = X * diag(1./std(X));

Xtest = Xtest - ones(size(Xtest,1),1)*mean(Xtest);
Xtest = Xtest * diag(1./std(Xtest));

% ... and plot
figure(3); clf
plot(X(Y==1,1),X(Y==1,2),'ro','MarkerSize',5);
hold on
plot(X(Y==-1,1),X(Y==-1,2),'b*','MarkerSize',5);
title('chessboard training data','FontSize',20);
axis([-2,2,-2,2]); 
axis square
xlabel('Dimension 1','FontSize',16);
ylabel('Dimension 2','FontSize',16);
set(gca,'FontSize',16);

% and wait for user to acknowledge that
fprintf('press key\n');
pause


%
% Start of MKL algorithm
%

% some initial guess of the optimal kernel width ...
dat = dist_euclid(X,X);
dat_test = dist_euclid(Xtest,X);
sig0 = 0.5*sqrt(median(dat(:)));

% ... and some widhts aroung this guess
sigs = sig0 * 10.^[-1:0.25:1];

% ... gram matrix for each width
for k=1:numel(sigs)
    Ks(:,:,k) = exp(-sigs(k)*dat);
    Ks_test(:,:,k) = exp(-sigs(k)*dat_test);
end
fprintf('using %d different kernels\n',size(Ks,3));

% train ...
C = 10;
[w,alpha,b] = simplemkl(Y,Ks,C);

% ... and test
Ytest_pred = sign(Kbeta(Ks_test,w) * alpha + b);


% print some statistics
fprintf('Algorithm selected %d kernels\n',sum(w>1e-4));
fprintf('test accuracy: %.5g\n',mean(Ytest==Ytest_pred));
[ignore,ind]= max(w);
fprintf('most active kernel (w=%.2g) = %.2f\n',w(ind),sigs(ind));



fprintf('press key to exit\n');
pause