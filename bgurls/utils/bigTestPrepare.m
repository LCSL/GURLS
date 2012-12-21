function [] = bigTestPrepare(filenameX, filenameY,files,blocksize)
% Prepare bigarrays for test set
% INPUT:
% -filenameX = file containg input data (nxd)
% -filenameY = file containing labels data (nxT for OVA nx1 otherwise)
% -files = structure containing the fields:
%     -Xtest_filename = prefix of files that make the bigarray for input
%       training data
%     -yest_filename = prefix of files that make the bigarray for ouput
%       training data
% -blocksize = number of samples per block


%% Test set.
Xte = bigarray_mat(files.Xtest_filename);
Xte.Clear();
Xte.Init(blocksize);

%% Test set labels.
yte = bigarray_mat(files.ytest_filename);
yte.Clear();
yte.Init(blocksize);

%% Read labels and assign samples to train, validation and test set balancing the classes
Yvec = importdata(filenameY);
n = length(Yvec);
Yvec = reshape(Yvec,n,1); %make sure it is a nx1 vector
T = max(Yvec); %number of classes

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
    
    Xte.Append(x);
    yte.Append(y);

    i_start = i_end +1;
    i_end = min(n,i_start + blocksize-1);
end	


fprintf('\n');

%% Make sure everything is written to disk.

Xte.Flush();
yte.Flush();

