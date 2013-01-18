% This demo uses the BIO data, stored in bio_unique.zip. 
% The data is in the form of an input data matrix and a outpit labels vector, 
% The data in split into training and test set.
% Then the bgurls training process is run on the training set and the bgurls testing
% process on the test set


%% Preprocess data
unzip('bio_unique.zip','bio_unique')
filenameX = 'bio_unique/X.csv'; %nxd input data matrix
filenameY = 'bio_unique/y.csv'; %nx1 or 1xn labels vector
blocksize = 1000; %matrices of size blocksizexd must fit into memory 
test_hoproportion = .2; %fraction of total samples to be used for testing
va_hoproportion = .2;  %fraction of training samples to be used for validation

dpath = 'bio_data_processed'; %direcory where all processed data is going to be stored

mkdir(dpath)

%set the prefix of the files that will constitute the bigarrays
files.Xtrain_filename = fullfile(dpath, 'bigarrays/Xtrain');
files.ytrain_filename = fullfile(dpath, 'bigarrays/ytrain');
files.Xtest_filename = fullfile(dpath, 'bigarrays/Xtest');
files.ytest_filename = fullfile(dpath, 'bigarrays/ytes');
files.Xva_filename = fullfile(dpath, 'bigarrays/Xva');
files.yva_filename = fullfile(dpath, 'bigarrays/yva');

%set the name of the files where pre-computed matrices will be stored
files.XtX_filename = fullfile(dpath, 'XtX.mat');
files.Xty_filename = fullfile(dpath, 'Xty.mat');
files.XvatXva_filename = fullfile(dpath,'XvatXva.mat');
files.Xvatyva_filename = fullfile(dpath, 'Xvatyva.mat');


% create bigarrays, by reading a unique input data file an a unique label data
% file and splitting the data into train, validation and test set.

fprintf('---preparing bigarrays...\n')
tic
bigTrainTestPrepare(filenameX, filenameY,files,blocksize,va_hoproportion,test_hoproportion)
toc
fprintf('---bigarrays prepared\n\n')

% compute and store matrices XtX, Xty, XvatXva, Xvatyva
tic
fprintf('---pre-computing relevant matrices...\n')
bigMatricesBuild(files)
toc
fprintf('---matrices computed \n\n')

%% Define the experiment options

name = 'bio_demoA';
opt = bigdefopt(name);

opt.files = files;

opt.files = rmfield(opt.files,{'Xtrain_filename';'ytrain_filename';'Xtest_filename';'ytest_filename'});
opt.files.pred_filename = fullfile(dpath, 'bigarrays/pred');


%% Define the learning pipeline

opt.seq = {'paramsel:dhoprimal','rls:dprimal','pred:primal','perf:macroavg'};
opt.process{1} = [2,2,0,0];
opt.process{2} = [3,3,2,2];

%% "Load" bigarrays variables for the training set
X = bigarray.Obj(files.Xtrain_filename);
y = bigarray.Obj(files.ytrain_filename);
X.Transpose(true);
y.Transpose(true);

%% Run bgurls on the training set
fprintf('---training...\n')
bgurls(X,y,opt,1)
fprintf('---training done\n\n')
%% "Load" bigarrays variables for the test set

X = bigarray.Obj(files.Xtest_filename);
y = bigarray.Obj(files.ytest_filename);
X.Transpose(true);
y.Transpose(true);

%% Run bgurls on the test set
fprintf('---testing...\n')
bgurls(X,y,opt,2);
fprintf('---testing done\n\n')

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
