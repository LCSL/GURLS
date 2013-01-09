function [] = bigTrainTestPrepare_manyfiles(DataDir,files,blocksize,test_hoproportion,val_hoproportion)
% Prepare bigarrays for training, validation and test set
% INPUT:
% -DataDir = directory cointaing as many .mat files as the number of classes 
%     Each file must contain the variables:
%     - x, the d x n_i input data matrix where d is the number of variables and n_i is the number of samples belonging to the i-th class;
% -files = structure containing the fields:
%     -Xtrain_filename = prefix of files that make the bigarray for input
%       training data
%     -ytrain_filename = prefix of files that make the bigarray for output
%       training data
%     -Xtest_filename = prefix of files that make the bigarray for input
%       test data
%     -ytest_filename = prefix of files that make the bigarray for output
%       test data
%     -Xva_filename = prefix of files that make the bigarray for input
%       validation data
%     -yva_filename = prefix of files that make the bigarray for output
%       validation data
% -blocksize = number of samples per block
% -test_hoproportion = percentage of total samples to be used for testing
% -val_hoproportion = percentage of training samples to be used for
% validation


list = dir(fullfile(DataDir, '*.mat'));
list = {list(:).name};
T = numel(list); %number of classes

%% Full training set.
Xtr = bigarray_mat(files.Xtrain_filename);
Xtr.Clear();
Xtr.Init(blocksize);

%% Full traning labels.
ytr = bigarray_mat(files.ytrain_filename);
ytr.Clear();
ytr.Init(blocksize);

%% Test set.
Xte = bigarray_mat(files.Xtest_filename);
Xte.Clear();
Xte.Init(blocksize);

%% Test set labels.
yte = bigarray_mat(files.ytest_filename);
yte.Clear();
yte.Init(blocksize);

%% Hold-out set (for cross validation), this is a subset of the training set.
Xva = bigarray_mat(files.Xva_filename);
Xva.Clear();
Xva.Init(blocksize);

%% Hold out set labels (subset ov ytr).
yva = bigarray_mat(files.yva_filename);
yva.Clear();
yva.Init(blocksize);

codes = 2*eye(T)-1;

%% Actually generate the dataset copying data in the bigarrays.

for i = 1:T
		load(fullfile(DataDir,list{i}));
        n_i = size(x,2);
		order = randperm(n_i);
		stop = floor(test_hoproportion*n_i);
		yTemp = repmat(codes(i,:),n_i,1);

		xte = x(:,order(1:stop));
		xtr = x(:,order(stop+1:end));
		

		Yte = yTemp(1:stop,:)';
		Ytr = yTemp(stop+1:end,:)';

		% take val_hoproportion of points for each class and throw them in the validation set.

		stop = floor(val_hoproportion*size(xtr,2));
		xva = xtr(:,1:stop);
		Yva = Ytr(:,1:stop);

		fprintf('Processing file %s: %.2f%%\r',list{i}, i*100/T);
		Xtr.Append(xtr);
		ytr.Append(Ytr);

		Xte.Append(xte);
		yte.Append(Yte);

		Xva.Append(xva);
		yva.Append(Yva);

end	
fprintf('\n');

%% Make sure everything is written to disk.

Xtr.Flush();
ytr.Flush();
Xte.Flush();
yte.Flush();
Xva.Flush();
yva.Flush();
