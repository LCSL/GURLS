% This demo uses the yeast data. 
% The data is already split into training and test set, and each set is 
% in the form of an input data matrix and a output labels vector. 
% Parameter selection and initial RLS estimation is carried out on a
% first subset of the training set. 
% Recursive RLS is run on the remainder of the training set, 
% simulating online learning.
% Finally the gurls testing process is run on the test set.

clear all
addpath('utils')
gurls_install
resdir = 'RESULTS';
mkdir(resdir);

n0 = 100; %size of first batch to be used for initialization
hoproportion = 0.1;

load(fullfile(gurls_root,'demo/data/yeast_data.mat'));
Xtr_tot = Xtr;
ytr_tot = ytr;
[ntr_tot,d] = size(Xtr_tot);


%% TRAIN: Use first batch for parameter selection and initialization

Xtr = Xtr_tot(1:n0,:);
ytr = ytr_tot(1:n0,:);

name = [resdir '/recursiveRLS'];
opt = gurls_defopt(name);
opt.kernel.XtX = Xtr'*Xtr;
opt.kernel.Xty = Xtr'*ytr;
opt.kernel.n = size(ytr,1);
opt.seq = {'split:ho','paramsel:hoprimal','rls:primalrecinit'};
opt.process{1} = [2,2,2];
opt = gurls(Xtr, ytr, opt,1);

Xva = Xtr(opt.split{1}.va,:);
yva = ytr(opt.split{1}.va,:);

%% UPDATE: Update RLS estimator recursively

for i = n0+1:ntr_tot    
    
    Xnew = Xtr_tot(i,:);
    ynew = ytr_tot(i,:);
    
    opt.rls = rls_primalrecupdate(Xnew, ynew, opt);
    
    opt.kernel.XtX = opt.kernel.XtX + Xnew'*Xnew;
    opt.kernel.Xty = opt.kernel.Xty + Xnew'*ynew;
    opt.kernel.n = opt.kernel.n + 1;

    if mod(opt.kernel.n,round(1/hoproportion))==0
        Xva = [Xva; Xnew];
        yva = [yva; ynew];
    end
end

%% EVAL: Test on independent test set
opt.pred = pred_primal(Xte, yte, opt); 
opt.perf = perf_macroavg(Xte, yte, opt); 

%% RETRAIN: retrain RLS and test

optRetrained = gurls_defopt(name);
optRetrained.kernel = opt.kernel;
nva = size(yva,1);
optRetrained.newprop('split', {});
optRetrained.split{1}.tr = zeros(opt.kernel.n - nva , 1);
optRetrained.split{1}.va = 1:nva;

optRetrained.seq = {'paramsel:hoprimal','rls:primalrecinit','pred:primal','perf:macroavg'};
optRetrained.process{1} = [2,2,0,0];
optRetrained.process{2} = [3,3,2,2];
optRetrained = gurls(Xva, yva, optRetrained,1);
optRetrained = gurls(Xte, yte, optRetrained,2);
 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('maximum difference between predicted Ys when using parameter chosen on the first batch or globally:')
disp(max(max(abs(opt.pred-optRetrained.pred))))

