clear all;
close all;

% dataset = 'breastcancer';
dataset = 'yeast';
if strcmp(dataset,'breastcancer')
    load(fullfile(gurls_root, 'demo/data/breastcancer_data.mat')); 
    ytr = ytr*2-1; yte = yte*2-1; % Need to change output to {-1,1}  
elseif strcmp(dataset,'yeast')
    load(fullfile(gurls_root, 'demo/data/yeast_data.mat'));
end

% filestr = {'linsvmsub','rbfsvmsub'};
filestr = {'linsvmsub','rbfsvmsub','rbfho','rbfloo'};
n = 5;

res_root = fullfile(gurls_root, 'demo/demodata'); % location where res files are stored

for r = 1:n
    strind = strmatch('linsvmsub',filestr);
    if ~isempty(strind)
        % Gaussian kernel, SVM Subgradient, Hold Out cross validation to select lambda and the Kernel width sigma
        name = [filestr{strind} '_' num2str(r)];
        opt = defopt(name);
        opt.seq = {'split:ho', 'kernel:linear', 'paramsel:hodual', 'rls:svmsubgradient', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
        opt.process{1} = [2,2,2,2,0,0,0];
        opt.process{2} = [3,3,3,3,2,2,2];
        opt.epochs = 100;
        gurls(Xtr, ytr, opt, 1);
        op = gurls(Xte, yte, opt, 2);
        disp(op.perf.acc);
    end

    strind = strmatch('rbfsvmsub',filestr);
    if ~isempty(strind)
        % Gaussian kernel, SVM Subgradient, Hold Out cross validation to select lambda and the Kernel width sigma
        name = [filestr{strind} '_' num2str(r)];
        opt = defopt(name);
        opt.seq = {'split:ho', 'paramsel:siglamho', 'kernel:rbf', 'rls:svmsubgradient', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
        opt.process{1} = [2,2,2,2,0,0,0,0];
        opt.process{2} = [3,3,3,3,2,2,2,2];
        opt.epochs = 100;
        gurls(Xtr, ytr, opt, 1);
        op = gurls(Xte, yte, opt, 2);
        disp(op.perf.acc);
    end

    strind = strmatch('rbfho',filestr);
    if ~isempty(strind)
        % Gaussian kernel, (dual formulation), Hold Out cross validation to select lambda and the Kernel width sigma
        name = [filestr{strind} '_' num2str(r)];
        opt = defopt(name);
        opt.seq = {'split:ho', 'paramsel:siglamho', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
        opt.process{1} = [2,2,2,2,0,0,0,0];
        opt.process{2} = [3,3,3,3,2,2,2,2];
        gurls(Xtr, ytr, opt, 1);
        gurls(Xte, yte, opt, 2); 
    end
    
    strind = strmatch('rbfloo',filestr);
    if ~isempty(strind)
        % Gaussian kernel, (dual formulation), Leave One Out cross validation to select lambda and the Kernel width sigma
        name = [filestr{strind} '_' num2str(r)];
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

nRuns = n*ones(1,length(filestr));
fields = {'perf.ap', 'perf.acc'};
plotopt.titles = {'Model Comparison - Accuracy', 'Model Comparison - Precision'};

if strcmp(dataset,'breastcancer')
    plotopt.class_names = {'Positive'};
elseif strcmp(dataset,'yeast')
    plotopt.class_names = {'CYT', 'NUC', 'MIT', 'ME3'};
end

% Generates "per-class" plots
summary_plot(filestr, fields, nRuns, plotopt, res_root)
% Generates "global" plots
summary_overall_plot(filestr, fields, nRuns, plotopt, res_root)
% Plots times taken by each step of the pipeline for performance reference.
plot_times(filestr, nRuns, res_root)

