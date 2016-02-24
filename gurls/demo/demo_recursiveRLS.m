% This demo uses the yeast data:  Data is already split into training and 
% test set, and each set is in the form of an input data matrix and a output 
% labels vector. 
% 
% - Parameter selection and initial RLS estimation is carried out on a
% first subset of the training set. 
%
% - Recursive RLS is run on the remainder of the training set, simulating 
% online learning.
%
% - GURLS testing process is run on the test set.

clear 
% addpath('utils')
% gurls_install
resdir = 'res_recursiveRLS';
if ~isdir(resdir), mkdir(resdir); end

load(fullfile(gurls_root,'demo/data/yeast_data.mat'));
Xtr_tot = Xtr;
ytr_tot = ytr;
[ntr_tot,d] = size(Xtr_tot);

n0 = 100; % size of first batch to be used for initialization

%% TRAIN: use first batch for parameter selection and initialization
Xtr = Xtr_tot(1:n0,:);
ytr = ytr_tot(1:n0,:);

name = [resdir '/recursiveRLS'];
opt = gurls_defopt(name);
opt.seq = {'split:ho','paramsel:hoprimal','rls:primalrecinit'};
opt.process{1} = [2,2,2];
opt = gurls(Xtr, ytr, opt,1);

%% UPDATE: update RLS estimator recursively
% update estimator recursively
for i = (n0+1):ntr_tot    
    
    Xnew = Xtr_tot(i,:);
    ynew = ytr_tot(i,:);
    
    opt.rls = rls_primalrecupdate(Xnew, ynew, opt);        
end

%% EVAL: test on independent test set
opt.pred = pred_primal(Xte, yte, opt); 
opt.perf = perf_macroavg(Xte, yte, opt); 

%% compare with batch RLS
name = [resdir '/standardRLSfixlambda'];
optRLSfixlambda = gurls_defopt(name);
optRLSfixlambda.newprop('paramsel.lambdas', opt.paramsel.lambdas*n0/ntr_tot);

optRLSfixlambda.seq = {'rls:primal','pred:primal','perf:macroavg'};
optRLSfixlambda.process{1} = [2,0,0];
optRLSfixlambda.process{2} = [3,2,2];
optRLSfixlambda = gurls(Xtr_tot, ytr_tot, optRLSfixlambda,1);
optRLSfixlambda = gurls(Xte, yte, optRLSfixlambda,2);

% check that solutions coincide
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('maximum difference between predicted Ys when using recursive and non recursive RLS:')
disp(max(max(abs(opt.pred-optRLSfixlambda.pred))))
