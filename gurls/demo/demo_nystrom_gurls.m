load(fullfile(gurls_root,'demo/data/breastcancer_data.mat'));

% Gaussian kernel approximated with Nystrom subsampling (m=200).
% Hold Out cross validation to select lambda and the Kernel width sigma

name = 'nysrbfho';
opt = gurls_defopt(name);
opt.seq = {'split:ho', 'paramsel:siglamho_nystrom', ...
    'kernel:rbf_nystrom', 'rls:dual_nystrom', ...
    'predkernel:traintest_nystrom', ...
    'pred:dual_nystrom', 'perf:rmse'};
opt.process{1} = [2,2,2,2,0,0,0];
opt.process{2} = [3,3,3,3,2,2,2];

opt.nystrom.shuffle = 0;
opt.nystrom.m = 200;
opt.hoperf = @perf_rmse;

gurls (Xtr, ytr, opt,1);
gurls (Xte, yte, opt,2);