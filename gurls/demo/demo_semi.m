%load(fullfile(gurls_root, 'demo/data/circles_data.mat'));
load(fullfile(gurls_root, 'demo/data/semi_data.mat'));
res_root = fullfile(gurls_root, 'demo'); % location where res files are stored

% Set up data (semi-supervised learning uses NaN for the label of unlabelled points)
nlabelled = 6; % number of labelled points
Xun = Xtr(nlabelled+1:1000,:);
Xtr = Xtr(1:nlabelled,:);  % use first few 'labelled' training points 
ytr = ytr(1:nlabelled,:);
Xtr2 = [Xtr; Xun]; % set the rest of the labels to NaN
ytr2 = [ytr; nan*ones(size(Xun,1),1)];

% Gaussian kernel, (dual formulation), Leave One Out cross validation to select lambda and the Kernel width sigma
name = 'rbfloocv';
opt = gurls_defopt(name);
opt.seq = {'paramsel:siglam', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg'};
opt.process{1} = [2,2,2,0,0,0];
opt.process{2} = [3,3,3,2,2,2];
opt.process{3} = [3,3,3,4,4,4];
opt = gurls(Xtr, ytr, opt, 1);
gurls(Xte, yte, opt, 2);

% Semi-supervised learning, Gaussian kernel, (dual formulation), Leave One Out cross validation
% (or fixed) to selection of lambdas and the Kernel widths
name = 'semiloocv';
opt2 = gurls_defopt(name);
%opt2.seq = {'paramsel:siglam_semi', 'kernel:rbf', 'rls:dual_semi', 'predkernel:traintest', 'pred:dual', 'perf:macroavg'};
opt2.seq = {'paramsel:fixsiglam_semi', 'kernel:rbf', 'rls:dual_semi', 'predkernel:traintest', 'pred:dual', 'perf:macroavg'};
opt2.process{1} = [2,2,2,0,0,0];
opt2.process{2} = [3,3,3,2,2,2];
opt2.process{3} = [3,3,3,4,4,4];
opt2.nlambda = 2;
opt2.nsigma = 2;
opt2 = gurls(Xtr2, ytr2, opt2, 1);
gurls(Xte, yte, opt2, 2);

%% Plot decision boundaries
v = -3:0.1:3;
ngrid = length(v);
[Xgr1,Xgr2]=meshgrid(v,v);
Xgr = [Xgr1(:),Xgr2(:)];
ygr = ones(length(Xgr),1);
XX = reshape(Xgr(:,1),ngrid,ngrid);
YY = reshape(Xgr(:,2),ngrid,ngrid);

% Gaussian kernel
fprintf('Gaussian kernel, accuracy:%8.4f\n',opt.perf.acc);
opt = gurls(Xgr, ygr, opt, 3);
opt = gurls(Xgr, ygr, opt, 2);
ZZ = reshape(opt.predkernel.K*opt.rls.C,ngrid,ngrid);

figure; 
[~,h] =contourf(XX,YY,ZZ);
axis equal; axis([-3 3 -3 3]); set(h,'edgecolor','none');
hold on;
scatter(Xun(:,1),Xun(:,2),'k.');
scatter(Xtr(ytr==1,1),Xtr(ytr==1,2),100,'rx','LineWidth',3);
scatter(Xtr(ytr==-1,1),Xtr(ytr==-1,2),100,'bx','LineWidth',3);
contour(XX,YY,ZZ,[0 0],'k','LineWidth',3);
cmap = jet;
colormap(0.2*cmap+0.8);
caxis([-1 1]);

% Semi-supervised learning
fprintf('Semi-supervised, accuracy:%8.4f\n',opt2.perf.acc);
opt2 = gurls(Xgr, ygr, opt2, 3);
opt2 = gurls(Xgr, ygr, opt2, 2);
ZZ = reshape(opt2.predkernel.K*opt2.rls.C,ngrid,ngrid);

figure;
[~,h] =contourf(XX,YY,ZZ);
axis equal; axis([-3 3 -3 3]); set(h,'edgecolor','none');
hold on;
scatter(Xun(:,1),Xun(:,2),'k.');
scatter(Xtr(ytr==1,1),Xtr(ytr==1,2),100,'rx','LineWidth',3);
scatter(Xtr(ytr==-1,1),Xtr(ytr==-1,2),100,'bx','LineWidth',3);
contour(XX,YY,ZZ,[0 0],'k','LineWidth',3);
colormap(0.2*cmap+0.8);
caxis([-1 1]);

