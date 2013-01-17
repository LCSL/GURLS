function [vout] = bigsplit_ho(X,y,opt)
% [splits] = bigsplit_ho(X, y, opt)
% Splits data into train and validation set

% INPUTS:
% -X: input data bigarray
% -y: labels bigarray
% -OPT: structure of options, with the following field with default values
%       set through the bigdefopt function:
%       -nholdouts
%       -hoproportion
% 
% OUTPUT: struct array of length OPT.NHOLDOUTS with each element gaving the following fields:
%   -tr: indices of samples to be used for training
%   -va: indices of samples to be used for testing

if isa(X,'bigarray')
	n = X.NumItems();
else	
	error('bigsplit:ho is to be used with bigarrays only');
end	

order = randperm(n);
last = floor(n*opt.hoproportion);
vout.va = order(1:last);
vout.tr = order(last+1:end);

%% Write the training bigarray X
%bXtr = bigarray_mat(opt.files.Xtr_filename);
%bXtr.Clear();
%bXtr.Init(X.BlockSize());
%X.Transpose(false);
%bXtr.Transpose(false);
%ba_csubset(X, bXtr, vout.tr);
%
%%% Write the training bigarray Y
%bYtr = bigarray_mat(opt.files.ytr_filename);
%bYtr.Clear();
%bYtr.Init(y.BlockSize());
%y.Transpose(false);
%bYtr.Transpose(false);
%ba_csubset(y, bYtr, vout.tr);


%% Write the validation bigarray

bXva = bigarray_mat(opt.files.Xva_filename);
bXva.Clear();
bXva.Init(X.BlockSize());
X.Transpose(false);
bXva.Transpose(false);
ba_csubset(X, bXva, vout.va);

bYva = bigarray_mat(opt.files.yva_filename);
bYva.Clear();
bYva.Init(y.BlockSize());
y.Transpose(false);
bYva.Transpose(false);
ba_csubset(y, bYva, vout.va);


% Reset the transpose.
X.Transpose(true);
y.Transpose(true);
