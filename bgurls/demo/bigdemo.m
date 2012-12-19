%=== WARNING===
% This script will download over 1 GB of data on your hard drive.

function [] = bigdemo()
    message = ['This script will donwload over 1 GB of data on your	hard drive, do you wish to continue? (y/[n]) : '];
    decision = userinput(message);

    if (~decision)
            return
    end

    %% Download data.
    fprintf('Downloading data... ')
    basePath = bgurls_root();
    dst = fullfile(basePath,'demo/data');
    untar('http://bratwurst.mit.edu/sbow.tar',dst);
    fprintf('Done!\n');


    %% Generate dataset

    % This converts ImageNet into bigarray using a 'standardized' directory
    % tree structure so that one can use bgTrainPrepare and bgTrainRun.
    %
    % If you wish to use a different directory tree structure just have a look at
    % bgTrainPrepare and bgTrainRun to see how gdm works.

    genDSet(dst, fullfile(dst,'ba'));

    %% Set up distributed matrix-matrix multiplications with gdm

    bgTrainPrepare(fullfile(dst,'ba'));

    %% Actually run distributed matrix-matrix multiplications with gdm (run this on multiple machines after you see the message!!).

    fprintf('Run this\n\tbgTrainRun(''%s'')\non as many machines as you please\n',fullfile(dst,'ba'));
    bgTrainRun(fullfile(dst,'ba'));

    %% "Load" bigarrays variables

    dpath = fullfile(dst,'ba/big');
    wpath = fullfile(dst,'ba');

    X = bigarray.Obj(fullfile(dpath,'trainX'));
    y = bigarray.Obj(fullfile(dpath,'trainY'));
    X.Transpose(true);
    y.Transpose(true);

    %% Define the experiment options

    % Names are pretty self-evident but just in case, have a look at bgTrainPrepare
    % to understand what these products are.

    name = fullfile(wpath,'gurls');
    opt = bigdefopt(name);
    opt.nb_pred = 5;

    opt.files.Xva_filename = fullfile(wpath, 'split/Xva');
    opt.files.yva_filename = fullfile(wpath, 'split/yva');

    opt.files.XtX_filename = fullfile(wpath, 'XtX.mat');
    opt.files.Xty_filename = fullfile(wpath, 'Xty.mat');


    opt.files.XvatXva_filename = fullfile(wpath,'XvatXva.mat');
    opt.files.Xvatyva_filename = fullfile(wpath, 'Xvatyva.mat');

    opt.files.pred_filename = fullfile(dpath, 'pred');

    %% Define the learning pipeline

    opt.seq = {'paramsel:dhoprimal','rls:dprimal','pred:primal','perf:macroavg'};
    opt.process{1} = [2,2,0,0];
    opt.process{2} = [3,3,2,2];

    %% Run bgurls on the training set

    bgurls(X,y,opt,1)

    X = bigarray.Obj(fullfile(dpath,'testX'));
    y = bigarray.Obj(fullfile(dpath,'testY'));
    X.Transpose(true);
    y.Transpose(true);

    %% Run bgurls on the test set

    bgurls(X,y,opt,2);

    % Now you should have a mat file in "wpath" named gurls.mat.
    % This file contains all the information about your experiment.
    % If you want to see the mean accuracy, for example, load the file
    % in your workspace and type
    %
    % >> mean(opt.perf.acc)
    %
    % If you are interested in visualizing or printing stats and facts
    % about your experiment, check the documentation about the summarizing
    % functions in the gurls package.
end

function [decision] = userinput(message)
choice = input(message,'s');
if isempty(choice) choice = 'n'; end
while ~strcmp(choice,'y') && ~strcmp(choice,'n')
    choice = input('\Please type ''y'' or ''n'' : ','s');
end
switch(choice)
    case {'y'}
        decision = true;
        return;
    case{'n'}
        decision = false;
        return;
end
end
