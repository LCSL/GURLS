%%%
%%% Example script for KPCA
%%%

%%
%
% Step 1: Make a data set.  This is similar to the 'circles' data set that
%  is used in the Wikipedia article for KPCA. I like it, so I will use it
%  here.
%
Ncirc = 5;
Nsamp = 100; % Per circle
% Make circles of data
X = zeros(0,2);


for k=1:Ncirc
    z = k^2*(0.5+0.2*rand(Nsamp,1)).*exp(2*pi*1i*rand(Nsamp,1));
    X = vertcat(X,[real(z) imag(z)]);
end
y = sqrt(sum(X.^2,2));

% (Linearlly) embed the data in a higher dimensional space.
U = orth(randn(50,2)); % Projector
X = X*U';

% Center & normalize
X = bsxfun(@rdivide,X,std(X));
y = y - mean(y);
y = y / std(y);


% Split data into testing and training sets
split = randperm(size(X,1));
M = floor(length(split)/2);
Xtr = X(split(1:M),:);
ytr = y(split(1:M));
Xte = X(split((M+1):end),:);
yte = y(split((M+1):end));

%% Unsupervised: Just do PCA and KPCA

% Setup options. This chain has no RLS in it.
opt = gurls_defopt('test');
opt.preproc.kernel.kernel = 'linear';
opt.preproc.n_dims = 2;
opt.seq = {'preproc:kpca_train','preproc:kpca_run'};
opt.process = {[2,2]};
[~,X_PCA] = gurls(X,[],opt,1);

opt.preproc.kernel.kernel = 'poly';
opt.preproc.kernel.c = 1;
opt.preproc.kernel.order = 2;
[~,X_KPCA] = gurls(X,[],opt,1);

figure(1);
subplot(1,2,1);
scatter(X_PCA(:,1),X_PCA(:,2),10,y);
axis equal;
title 'Data after PCA';
subplot(1,2,2);
scatter(X_KPCA(:,1),X_KPCA(:,2),10,y);
axis equal;
title 'Data after KPCA (poly kernel)';
%% Supervised: Regression model w/ dimensionality selection via crossvalidation
opt = gurls_defopt('test');
opt.preproc.kernel.kernel = 'poly';
opt.preproc.kernel.c = 1;
opt.preproc.kernel.order = 2;
opt.preproc.n_dims = 232; % set to whatever, this is set by cross-validation.
opt.kernel.type = 'linear';
opt.verbose = false;
opt.seq = {'preproc:kpca_train','preproc:kpca_run','rls:primal','pred:primal','perf:rmse'};
opt.process = {[2,2,2,0,0],[3,2,2,2,2]};
opt.newprop('paramsel',struct());
opt.paramsel.lambdas = 1E-5;
cv_params = struct;
cv_params.n_trials = 10;
cv_params.pct_train = 0.75;
cv_params.param_name = 'preproc.n_dims';
cv_params.param_vals = 1:5;%[0 1 2 3 4 5];
cv_params.perf_field = 'rmse';
cv_params.max_or_min = 'min';

preproc_crossvalidate(Xtr,Xte,opt,cv_params);

gurls(Xtr, ytr, opt, 1);
[~,X_out]=gurls(Xte, yte, opt, 2);

figure(2); clf;
scatter(yte,opt.pred);
xlabel 'Truth value';
ylabel 'Model prediction';
title 'Linear regression after PCA';

