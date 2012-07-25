if ~ischar(datadir)
    error('datadir is not a string')
end
if isempty(datadir)
    error('datadir is empty')
end

cd(datadir)

%% default options
fid = fopen('singlelambda.txt','w');
fprintf(fid,'median');
fclose(fid);
dlmwrite('nlambda.txt',10);
dlmwrite('nsigma.txt',5);
dlmwrite('nholdouts.txt',1);
dlmwrite('hoproportion.txt',0.2);
dlmwrite('smallnumber.txt',1e-8,'precision','%1.11e');
fid = fopen('hoperf.txt','w');
fprintf(fid,'macroavg');
fclose(fid);

%% split
dirname = 'splitho';
if ~exist(dirname,'dir'); mkdir(dirname); end
test_driver('split_ho',dirname)

%% paramsel A
dirname = 'paramselsiglam';
if ~exist(dirname,'dir'); mkdir(dirname); end
test_driver('paramsel_siglam',dirname)

dirname = 'paramselsiglamho';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('splitho/split*',dirname)
test_driver('paramsel_siglamho',dirname)

%% kernel
dirname = 'kernellinear';
if ~exist(dirname,'dir'); mkdir(dirname); end
test_driver('kernel_linear',dirname)

dirname = 'kernelgauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('paramselsiglam/paramsel-sigma.txt',dirname)
test_driver('kernel_rbf',dirname)

dirname = 'kernelchisquared';
if ~exist(dirname,'dir'); mkdir(dirname); end
test_driver('kernel_chisquared',dirname)


%% paramsel B
dirname = 'paramselloocvprimal';
if ~exist(dirname,'dir'); mkdir(dirname); end
test_driver('paramsel_loocvprimal',dirname)

dirname = 'paramselhoprimal';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('splitho/split*',dirname)
test_driver('paramsel_hoprimal',dirname)

dirname = 'paramselloocvdual_linear';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('kernellinear/kernel-K.txt',dirname)
copyfile('kernellinear/kernel-type.txt',dirname)
test_driver('paramsel_loocvdual',dirname)

dirname = 'paramselhodual_linear';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('splitho/split*',dirname)
copyfile('kernellinear/kernel-K.txt',dirname)
copyfile('kernellinear/kernel-type.txt',dirname)
test_driver('paramsel_hodual',dirname)

dirname = 'paramselloocvdual_gauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('kernelgauss/kernel-K.txt',dirname)
copyfile('kernelgauss/kernel-type.txt',dirname)
test_driver('paramsel_loocvdual',dirname)

dirname = 'paramselhodual_gauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('splitho/split*',dirname)
copyfile('kernelgauss/kernel-K.txt',dirname)
copyfile('kernelgauss/kernel-type.txt',dirname)
test_driver('paramsel_hodual',dirname)

dirname = 'paramselhoprimalr';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('splitho/split*',dirname)
test_driver('paramsel_hoprimalr',dirname)

dirname = 'paramselhodualr_linear';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('splitho/split*',dirname)
copyfile('kernellinear/kernel-K.txt',dirname)
copyfile('kernellinear/kernel-type.txt',dirname)
test_driver('paramsel_hodualr',dirname)

dirname = 'paramselhodualr_gauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('splitho/split*',dirname)
copyfile('kernelgauss/kernel-K.txt',dirname)
copyfile('kernelgauss/kernel-type.txt',dirname)
test_driver('paramsel_hodualr',dirname)

%% rls

dirname = 'rlsprimal';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('paramselhoprimal/paramsel-lambdas.txt',dirname)
test_driver('rls_primal',dirname)

dirname = 'rlsdual_linear';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('paramselhodual_linear/paramsel-lambdas.txt',dirname)
copyfile('kernellinear/kernel-K.txt',dirname)
copyfile('kernellinear/kernel-type.txt',dirname)
test_driver('rls_dual',dirname)

dirname = 'rlsdual_gauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('paramselhodual_gauss/paramsel-lambdas.txt',dirname)
copyfile('kernelgauss/kernel-K.txt',dirname)
copyfile('kernelgauss/kernel-type.txt',dirname)
test_driver('rls_dual',dirname)

dirname = 'rlsprimalr';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('paramselhoprimalr/paramsel-lambdas.txt',dirname)
test_driver('rls_primalr',dirname)

dirname = 'rlsdualr_linear';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('paramselhodualr_linear/paramsel-lambdas.txt',dirname)
copyfile('kernellinear/kernel-K.txt',dirname)
copyfile('kernellinear/kernel-type.txt',dirname)
test_driver('rls_dualr',dirname)

dirname = 'rlsdualr_gauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('paramselhodualr_gauss/paramsel-lambdas.txt',dirname)
copyfile('kernelgauss/kernel-K.txt',dirname)
copyfile('kernelgauss/kernel-type.txt',dirname)
test_driver('rls_dualr',dirname)

%% predkernel

dirname = 'predkernelchisquared';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('kernelchisquared/kernel-type.txt',dirname)
copyfile('rlsdual_gauss/optimizer-X.txt',dirname) %va bene anche di un altro kernel, tanto prendo solo X
test_driver('predkernel_traintest',dirname)

dirname = 'predkernel_gauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('kernelgauss/kernel-type.txt',dirname)
copyfile('rlsdual_gauss/optimizer-X.txt',dirname) %va bene anche di un altro kernel, tanto prendo solo X
copyfile('paramselsiglam/paramsel-sigma.txt',dirname)
test_driver('predkernel_traintest',dirname)

%% pred

dirname = 'predprimal';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('rlsprimal/optimizer*',dirname)
test_driver('pred_primal',dirname)

dirname = 'preddual_linear';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('rlsdual_linear/optimizer*',dirname)
copyfile('kernellinear/kernel-type.txt',dirname)
test_driver('pred_dual',dirname)

dirname = 'preddual_gauss';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('rlsdual_gauss/optimizer*',dirname)
copyfile('predkernel_gauss/predkernel-K.txt',dirname)
copyfile('kernelgauss/kernel-type.txt',dirname)
test_driver('pred_dual',dirname)

%% perf

dirname = 'perfrmse';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('predprimal/pred*',dirname)
test_driver('perf_rmse',dirname)

dirname = 'perfmacroavg';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('predprimal/pred*',dirname)
test_driver('perf_macroavg',dirname)

dirname = 'perfprecrec';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('predprimal/pred*',dirname)
test_driver('perf_precrec',dirname)


%% conf

dirname = 'confgap';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('predprimal/pred*',dirname)
test_driver('conf_gap',dirname)

dirname = 'confmaxscore';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('predprimal/pred*',dirname)
test_driver('conf_maxscore',dirname)

dirname = 'confboltzman';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('predprimal/pred*',dirname)
test_driver('conf_boltzman',dirname)

dirname = 'confboltzmangap';
if ~exist(dirname,'dir'); mkdir(dirname); end
copyfile('predprimal/pred*',dirname)
test_driver('conf_boltzmangap',dirname)


