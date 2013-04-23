% This demo uses the yeast data. 
% The data is already split into training and test set, and each set is 
% in the form of an input data matrix and a output labels vector. 
% Parameter selection and initial RLS estimation is carried out on a
% first subset of the training set. 
% Recursive RLS is run on the remainder of the training set, 
% simulating online learning.
% Finally the gurls testing process is run on the test set.

clear all
resdir = 'RESULTS';
addpath('utilities')
gurls_install

%% load data

n0 = 100; %size of first batch to be used for initialization
epochs = 10;  % number of samples for upadte (a corresponding number of 
                % rows is read from the file where the global training set 
                % is stored (can also be epochs=1)
nval = 100; % number of samples to be used for validation in retraining

load(fullfile(gurls_root,'demo/data/yeast_data.mat'));
Xtr_tot = Xtr;
ytr_tot = ytr;
[ntr_tot,d] = size(Xtr_tot);


%% load first batch for parameter selection and initialization

mkdir(resdir);

Xtr = Xtr_tot(1:n0,:);
ytr = ytr_tot(1:n0,:);

name = [resdir '/recursiveRLS'];
opt = defopt(name);
opt.kernel.XtX = Xtr'*Xtr;
opt.kernel.Xty = Xtr'*ytr;
opt.kernel.n = size(ytr,1);
opt.seq = {'split:ho','paramsel:hoprimal','rls:primalrecinit'};
opt.process{1} = [2,2,2];
opt = gurls (Xtr, ytr, opt,1);

%% update RLS estimator recursively
% read a block from the file where the global training set is stored and
% update estimator
i_start = n0+1;
i_end = n0+epochs;
while i_start<ntr_tot-1    
    
    Xnew = Xtr_tot(i_start:i_end,:);
    ynew = ytr_tot(i_start:i_end,:);
    
    opt.rls = rls_primalrecupdate(Xnew, ynew, opt);
    
    opt.kernel.XtX = opt.kernel.XtX + Xnew'*Xnew;
    opt.kernel.Xty = opt.kernel.Xty + Xnew'*ynew;
    opt.kernel.n = opt.kernel.n + size(ynew,1);
    
    i_start = i_end+1;
    i_end = min(ntr_tot,i_end+epochs);
end

%% test on independent test set
opt.pred = pred_primal(Xte, yte, opt); 
opt.perf = perf_macroavg(Xte, yte, opt); 

%% retrain RLS and test

optRetrained = defopt(name);
optRetrained.kernel = opt.kernel;

val_indices = randperm(ntr_tot);
val_indices = val_indices(1:nval);
Xva = Xtr_tot(val_indices,:);
yva = ytr_tot(val_indices,:);
optRetrained.split.tr = zeros(opt.kernel.n,1);
optRetrained.split.va = 1:size(yva,1);

optRetrained.seq = {'paramsel:hoprimal','rls:primalrecinit','pred:primal','perf:macroavg'};
optRetrained.process{1} = [2,2,0,0];
optRetrained.process{2} = [3,3,2,2];
optRetrained = gurls (Xva, yva, optRetrained,1);
optRetrained = gurls (Xte, yte, optRetrained,2);


%%
name = [resdir '/standardRLSfixlambda'];
optRLSfixlambda = defopt(name);
optRLSfixlambda.paramsel.lambdas = opt.paramsel.lambdas;
optRLSfixlambda.seq = {'rls:primal','pred:primal','perf:macroavg'};
optRLSfixlambda.process{1} = [2,0,0];
optRLSfixlambda.process{2} = [3,2,2];
optRLSfixlambda = gurls (Xtr_tot, ytr_tot, optRLSfixlambda,1);
optRLSfixlambda = gurls (Xte, yte, optRLSfixlambda,2);

% check that solutions coincide
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('maximum difference between predicted Ys when using recursive and non recursive RLS:')
disp(max(max(abs(opt.pred-optRLSfixlambda.pred))))
disp(mean(opt.perf.acc)-mean(optRLSfixlambda.perf.acc))

%% what changes if the regularization is chosen globally or on the first batch?
name = [resdir '/standardRLS'];
optRLS = defopt(name);
optRLS.seq = {'split:ho','paramsel:hoprimal','rls:primal','pred:primal','perf:macroavg'};
optRLS.process{1} = [2,2,2,0,0];
optRLS.process{2} = [0,0,3,2,2];
optRLS = gurls (Xtr_tot, ytr_tot, optRLS,1);
optRLS = gurls (Xte, yte, optRLS,2);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('maximum difference between predicted Ys when using recursive and non recursive RLS:')
disp(max(max(abs(opt.pred-optRLS.pred))))
disp(mean(opt.perf.acc)-mean(optRLS.perf.acc))