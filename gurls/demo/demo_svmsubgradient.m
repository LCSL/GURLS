% DEMO OF THE SVM SUBGRADIENT METHOD COMPARED TO DUAL FORMULATION
% adapted from gurlsdemo.m

% Use either the 'breastcancer' or 'yeast' datasets
dataset = 'yeast';
if strcmp(dataset,'breastcancer')
    load(fullfile(gurls_root, 'demo/data/breastcancer_data.mat')); 
    ytr = ytr*2-1; yte = yte*2-1; % Convert output to {-1,1}  
elseif strcmp(dataset,'yeast')
    load(fullfile(gurls_root, 'demo/data/yeast_data.mat'));
end

% List the methods to use
models = {'linsvmsub','rbfsvmsub','rbfho','rbfloo'};

% Number of trials to run each method
n = 5;

res_root = fullfile(gurls_root, 'demo'); % location where res files are stored

for r = 1:n
    strind = strmatch('linsvmsub',models);
    if ~isempty(strind)
        % Gaussian kernel, SVM Subgradient, Hold Out cross validation to select lambda
        name = [models{strind} '_' num2str(r)];
        opt = defopt(name);
        opt.seq = {'split:ho', 'kernel:linear', 'paramsel:hodual', 'rls:svmsubgradient', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
        opt.process{1} = [2,2,2,2,0,0,0];
        opt.process{2} = [3,3,3,3,2,2,2];
        opt.Niter = 200;
%         opt.gammafunc = @(x) power(x,-.5);
        gurls(Xtr, ytr, opt, 1);
        gurls(Xte, yte, opt, 2);
    end

    strind = strmatch('rbfsvmsub',models);
    if ~isempty(strind)
        % Gaussian kernel, SVM Subgradient, Hold Out cross validation to select lambda and the Kernel width sigma
        name = [models{strind} '_' num2str(r)];
        opt = defopt(name);
        opt.seq = {'split:ho', 'paramsel:siglamho', 'kernel:rbf', 'rls:svmsubgradient', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
        opt.process{1} = [2,2,2,2,0,0,0,0];
        opt.process{2} = [3,3,3,3,2,2,2,2];
        opt.Niter = 200;
%         opt.gammafunc = @(x) power(x,-.5);
        gurls(Xtr, ytr, opt, 1);
        gurls(Xte, yte, opt, 2);
    end

    strind = strmatch('rbfho',models);
    if ~isempty(strind)
        % Gaussian kernel, (dual formulation), Hold Out cross validation to select lambda and the Kernel width sigma
        name = [models{strind} '_' num2str(r)];
        opt = defopt(name);
        opt.seq = {'split:ho', 'paramsel:siglamho', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
        opt.process{1} = [2,2,2,2,0,0,0,0];
        opt.process{2} = [3,3,3,3,2,2,2,2];
        gurls(Xtr, ytr, opt, 1);
        gurls(Xte, yte, opt, 2); 
    end
    
    strind = strmatch('rbfloo',models);
    if ~isempty(strind)
        % Gaussian kernel, (dual formulation), Leave One Out cross validation to select lambda and the Kernel width sigma
        name = [models{strind} '_' num2str(r)];
        opt = defopt(name);
        opt.seq = {'paramsel:siglam', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
        opt.process{1} = [2,2,2,0,0,0,0];
        opt.process{2} = [3,3,3,2,2,2,2];
        gurls(Xtr, ytr, opt, 1);
        gurls(Xte, yte, opt, 2);
    end
end

%% Visualization:
%	- filestr: cell with names of the experiments
%	- nRuns: number of runs for each experiment
%	- fields: which fields of opt to display (as many plots as the elements of fields will be generated).
%	- plotopt: a structure containing various text labels for the plots.

nRuns = n*ones(1,size(models,2));
fields = {'perf.ap', 'perf.acc'};
plotopt.titles = {'Model Comparison - Accuracy', 'Model Comparison - Precision'};

% Label based on dataset used
if strcmp(dataset,'breastcancer')
    plotopt.class_names = {'Positive'};
elseif strcmp(dataset,'yeast')
    plotopt.class_names = {'CYT', 'NUC', 'MIT', 'ME3'};
end

% Generates "per-class" plots
summary_plot(models, fields, nRuns, plotopt, res_root)
% Generates "global" plots
summary_overall_plot(models, fields, nRuns, plotopt, res_root)
% Plots times taken by each step of the pipeline for performance reference.
plot_times(models, nRuns, res_root)

