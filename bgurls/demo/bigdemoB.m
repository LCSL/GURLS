% This demo uses the BIO data, stored in bio_TrainTest.zip. 
% The data is in the form of two input data matrices and two outpit labels vectors,
% one for traning and one for test.
% Then run the bgurls training process on the training set and the bgurls testing
% process on the test set

%% Preprocess data
unzip('bio_TrainTest.zip','bio_TrainTest')
filenameXtrain = 'bio_TrainTest/Xtr.csv'; %nxd input data matrix for training
filenameYtrain = 'bio_TrainTest/ytr.csv'; %nx1 or 1xn labels vector for training
filenameXtest = 'bio_TrainTest/Xte.csv'; %nxd input data matrix for test
filenameYtest = 'bio_TrainTest/yte.csv'; %nx1 or 1xn labels vector for test
blocksize = 1000; %matrices of size blocksizexd must fit into memory 
test_hoproportion = .2; %fraction of total samples to be used for testing
va_hoproportion = .2;  %fraction of training samples to be used for validation

dpath = 'bio_data_processed'; %direcory where all processed data is going to be stored

mkdir(dpath)

%set the prefix of the files that will constitute the bigarrays (each
%bigarray is made of a set of file with the same prefix)
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


% create bigarrays for training and validation, by reading a unique input data file an a unique label data
% file and splitting the data into train and validation set.
fprintf('---preparing bigarrays for traing...\n')
tic
bigTrainPrepare(filenameXtrain, filenameYtrain,files,blocksize,va_hoproportion)
toc
fprintf('---training bigarrays prepared\n\n')

% create bigarrays for test
fprintf('---preparing bigarrays for test...\n')
tic
bigTestPrepare(filenameXtest, filenameYtest,files,blocksize)
toc
fprintf('---test bigarrays prepared\n\n')

% compute and store matrices XtX, Xty, XvatXva, Xvatyva
tic
fprintf('---pre-computing relevant matrices...\n')
bigMatricesBuild(files)
toc
fprintf('---matrices computed \n\n')

%% Define the experiment options

name = 'bio_demoB';
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
