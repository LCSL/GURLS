function [vout] = bigsplit_ho(X,y,opt)

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
