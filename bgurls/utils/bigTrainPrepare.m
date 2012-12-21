function [] = bigTrainPrepare(filenameX, filenameY,files,blocksize,val_hoproportion)
% Prepare bigarrays for training and validation
% INPUT:
% -filenameX = file containg input data (nxd)
% -filenameY = file containing labels data (nxT for OVA nx1 otherwise)
% -files = structure containing the fields:
%     -Xtrain_filename = prefix of files that make the bigarray for input
%       training data
%     -ytrain_filename = prefix of files that make the bigarray for ouput
%       training data
%     -Xva_filename = prefix of files that make the bigarray for input
%       validation data
%     -ytva_filename = prefix of files that make the bigarray for ouput
%       validation data
% -blocksize = number of samples per block
% -val_hoproportion = percentage of training samples to be used for validation


%% Full training set.
Xtr = bigarray_mat(files.Xtrain_filename);
Xtr.Clear();
Xtr.Init(blocksize);

%% Full traning labels.
ytr = bigarray_mat(files.ytrain_filename);
ytr.Clear();
ytr.Init(blocksize);

%% Hold-out set (for cross validation), this is a subset of the training set.
Xva = bigarray_mat(files.Xva_filename);
Xva.Clear();
Xva.Init(blocksize);

%% Hold out set labels (subset ov ytr).
yva = bigarray_mat(files.yva_filename);
yva.Clear();
yva.Init(blocksize);

%% Read labels and assign samples to train, validation and test set balancing the classes
Yvec = importdata(filenameY);
n = length(Yvec);
Yvec = reshape(Yvec,n,1); %make sure it is a nx1 vector

T = max(Yvec); %number of classes
val_indices = zeros(n,1);
for t = 1:T;
    class_indices = find(Yvec==t); 
    class_size = length(class_indices);
    perm = randperm(class_size);
    val_end = round(class_size*(val_hoproportion));
    val_indices_tmp = class_indices(perm(1:val_end));
    val_indices(val_indices_tmp) = 1;
        
end
val_indices = logical(val_indices);

%% Actually generate the dataset copying data in the bigarrays.
x = csvread(filenameX,n-1,0);
d = length(x);

i_start = 1;
i_end = blocksize;
codes = 2*eye(T) - 1;
while i_start<=n;
    % read blocksize rows of filenameX
    range = [i_start-1,0,i_end-1,d-1];
    x = csvread(filenameX,i_start-1,0,range)';
    
    block_indices = i_start:i_end;

    y = codes(Yvec(block_indices),:)';
    
    Xtr.Append(x);
    ytr.Append(y);

    % take 20% of points for each class and throw them in the validation set.
    Xva.Append(x(:,val_indices(block_indices)));
    yva.Append(y(:,val_indices(block_indices)));

    i_start = i_end +1;
    i_end = min(n,i_start + blocksize-1);
end	


fprintf('\n');

%% Make sure everything is written to disk.

Xtr.Flush();
ytr.Flush();
Xva.Flush();
yva.Flush();

