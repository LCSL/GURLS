clear; close all;
load(fullfile(gurls_root, 'demo/data/yeast_data.mat'));
% load(fullfile(gurls_root, 'demo/data/semi_data.mat'));
res_root = fullfile(gurls_root, 'demo'); % location where res files are stored

nTrials = 5;
%% NOTES:
% - Executes "nTrials" times 5 different pipelines on the same dataset.
% - Uses the summary routines to visualize precision and accuracy per class
% and globally for each pipeline.

for r = 1:nTrials
    
    % Linear kernel, primal formulation, Leave One Out cross validation to select lambda
    name = fullfile(res_root, ['linloo_' num2str(r)]);
    opt = gurls_defopt(name);
    opt.seq = {'paramsel:loocvprimal','rls:primal','pred:primal', 'perf:macroavg', 'perf:precrec'};
    opt.process{1} = [2,2,0,0,0];
    opt.process{2} = [3,3,2,2,2];
    gurls(Xtr, ytr, opt, 1);
    gurls(Xte, yte, opt, 2);
    
    % plot_decision_boundaries(Xtr, ytr, opt);
    
    % Linear kernel, dual formulation, Leave One Out cross validation to select lambda
    name = fullfile(res_root, ['lindloo_' num2str(r)]);
    opt = gurls_defopt(name);
    opt.seq = {'kernel:linear', 'paramsel:loocvdual', 'rls:dual', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
    opt.process{1} = [2,2,2,0,0,0];
    opt.process{2} = [3,3,3,2,2,2];
    gurls(Xtr, ytr, opt, 1);
    gurls(Xte, yte, opt, 2);
    
    % plot_decision_boundaries(Xtr, ytr, opt);
    
    % Gaussian kernel, Leave One Out cross validation to select lambda and the Kernel width sigma
    name = fullfile(res_root, ['rbfloo_' num2str(r)]);
    opt = gurls_defopt(name);
    opt.seq = {'paramsel:siglam', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
    opt.process{1} = [2,2,2,0,0,0,0];
    opt.process{2} = [3,3,3,2,2,2,2];
    % opt.nsigma = 10;
    % opt.nlambda = 10;
    gurls(Xtr, ytr, opt, 1);
    gurls(Xte, yte, opt, 2);
        
    % plot_decision_boundaries(Xtr, ytr, opt);
    
    % Linear kernel, primal formulation, Hold Out cross validation to select lambda
    name = fullfile(res_root, ['linho_' num2str(r)]);
    opt = gurls_defopt(name);
    opt.seq = {'split:ho','paramsel:hoprimal','rls:primal','pred:primal','perf:macroavg','perf:precrec'};
    opt.process{1} = [2,2,2,0,0,0];
    opt.process{2} = [3,3,3,2,2,2];
    gurls(Xtr, ytr, opt, 1);
    gurls(Xte, yte, opt, 2);
    
    % plot_decision_boundaries(Xtr, ytr, opt);
    
    % Gaussian kernel, Hold Out cross validation to select lambda and the Kernel width sigma
    name = fullfile(res_root, ['rbfho_' num2str(r)]);
    opt = gurls_defopt(name);
    opt.seq = {'split:ho', 'paramsel:siglamho', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
    opt.process{1} = [2,2,2,2,0,0,0,0];
    opt.process{2} = [3,3,3,3,2,2,2,2];
    gurls(Xtr, ytr, opt, 1);
    gurls(Xte, yte, opt, 2);
    
    % plot_decision_boundaries(Xtr, ytr, opt);
    
end

%% Visualization:
%	- filestr: cell with names of the experiments
%	- nRuns: number of runs for each experiment
%	- fields: which fields of opt to display (as many plots as the elements of fields will be generated).
%	- plotopt: a structure containing various text labels for the plots.

filestr = {'linloo', 'lindloo', 'rbfloo', 'linho', 'rbfho'};
nRuns = nTrials*ones(1,5);
fields = {'perf.ap', 'perf.acc'};
plotopt.titles = {'Model Comparison - Accuracy', 'Model Comparison - Precision'};
plotopt.class_names = {'CYT', 'NUC', 'MIT', 'ME3'};

% Generates "per-class" plots
summary_plot(filestr, fields, nRuns, plotopt, res_root)
% Generates "global" plots
summary_overall_plot(filestr, fields, nRuns, plotopt, res_root)
% Plots times taken by each step of the pipeline for performance reference.
plot_times(filestr, nRuns, res_root)

