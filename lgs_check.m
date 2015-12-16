%basic example: linear kernel
%simulate data
X = normrnd(0,1,1000,2);
beta = [2;-1];
p = 1./(1+exp(-X*beta));
y = binornd( 1, p );
y(y==0) = -1;

name = 'BasicL'
opt = defopt(name);
opt.seq = ...
{'split:ho','paramsel:bfprimal','rls:landweberprimal_lgs', ...
 'pred:primal','perf:macroavg'};
%{'split:ho','kernel:rbf','paramsel:bfdual','rls:landweberdual', ...
% 'predkernel:traintest','pred:dual','perf:macroavg'};
opt.process{1} = [2,2,2,0,0];
opt.process{2} = [3,3,3,2,2];
%For selecting optimal steps
opt.nholdouts = 5;
opt.paramsel.guesses = linspace( 100, 1000, 301 );
opt.paramsel.optimizer = @rls_landweberprimal_lgs;

res1 = gurls (X, y, opt, 1);
display( res1.rls.W );

%basic example: Gaussina kernel
%simulate data
X = normrnd(0,1,5000,2);
p = 1./(1+exp(-(X(:,2)-sin(X(:,1)))));
y = binornd( 1, p );
y(y==0) = -1;

name = 'BasicG'
opt = defopt(name);
opt.seq = ...
{'split:ho','kernel:rbf','paramsel:bfdual','rls:landweberdual', ...
 'predkernel:traintest','pred:dual','perf:macroavg'};
opt.process{1} = [2,2,2,2,0,0,0];
%For selecting optimal steps
opt.nholdouts = 5;
opt.paramsel.guesses = linspace( 100, 1000, 301 );
opt.paramsel.optimizer = @rls_landweberdual_lgs;
opt.paramsel.sigma = 1;

res2 = gurls (X, y, opt, 1);
res_pred = res2.kernel.K*res2.rls.C;

%plot the decision boundary
xlin = linspace( min(X(:,1)),max(X(:,1)),33 );
ylin = linspace( min(X(:,2)),max(X(:,2)),33 );
[Xmesh,Ymesh] = meshgrid( xlin, ylin );
f = scatteredInterpolant( X(:,1), X(:,2), res_pred );
Z = f(Xmesh,Ymesh);
[C,h] = contour( Xmesh, Ymesh, Z, [0,0] );
h.LineColor = 'k'
hold on
plot( linspace(min(X(:,1)),max(X(:,1)),100), sin(linspace(min(X(:,1)),max(X(:,1)),100)), 'r' );
plot( X(y==1,1), X(y==1,2), '.b', 'markers', 6 );
plot( X(y==-1,1), X(y==-1,2), '.g', 'markers', 6 );
legend('Estimated boundary','True boundary','Class 1','Class -1');
xlim([-3 3]);
ylim([-3 3]);
hold off

% Gaussian kernel, homework data
load('C:/Users/Renboyu/Desktop/9520_fall2015_pset1/9520_fall2015_pset1/2015_ps1-dataset_1.mat')
name = 'HomeworkG'
opt = defopt(name);
opt.seq = ...
{'split:ho','kernel:rbf','paramsel:bfdual','rls:landweberdual', ...
 'predkernel:traintest','pred:dual','perf:macroavg'};
opt.process{1} = [2,2,2,2,0,0,0];
opt.process{2} = [3,3,3,3,2,2,2];

%For selecting optimal steps
opt.nholdouts = 5;
opt.paramsel.guesses = linspace( 100, 1000, 301 );
opt.paramsel.optimizer = @rls_landweberdual_lgs;
opt.paramsel.sigma = mean( pdist( Xtrain ) );

gurls (Xtrain, Ytrain, opt, 1);
res3 = gurls( Xtest, Ytest, opt, 2);
display( res3.perf.acc );

%plot the decision boundary
xlin = linspace( min(Xtest(:,1)),max(Xtest(:,1)),33 );
ylin = linspace( min(Xtrain(:,2)),max(Xtrain(:,2)),33 );
[Xmesh,Ymesh] = meshgrid( xlin, ylin );
f = scatteredInterpolant( Xtest(:,1), Xtest(:,2), res3.pred );
Z = f(Xmesh,Ymesh);
[C,h] = contour( Xmesh, Ymesh, Z, [0,0] );
h.LineColor = 'k'
hold on
plot( Xtest(Ytest==1,1), Xtest(Ytest==1,2), '.b', 'markers', 6 );
plot( Xtest(Ytest==-1,1), Xtest(Ytest==-1,2), '.g', 'markers', 6 );
legend('Estimated boundary','Class 1','Class -1');
hold off