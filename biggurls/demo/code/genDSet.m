function dst = genDSet(src,dst)

list = dir(fullfile(src, 'n*.mat'));
list = {list(:).name};

wpath = dst;

%% Full training set.
Xtr = bigarray_mat(fullfile(wpath,'big/trainX'));
Xtr.Clear();
Xtr.Init(1000);

%% Full traning labels.
ytr = bigarray_mat(fullfile(wpath,'big/trainY'));
ytr.Clear();
ytr.Init(1000);

%% Test set.
Xte = bigarray_mat(fullfile(wpath,'/big/testX'));
Xte.Clear();
Xte.Init(1000);

%% Test set labels.
yte = bigarray_mat(fullfile(wpath,'/big/testY'));
yte.Clear();
yte.Init(1000);

%% Hold-out set (for cross validation), this is a subset of the training set.
Xva = bigarray_mat(fullfile(wpath,'split/Xva'));
Xva.Clear();
Xva.Init(1000);

%% Hold out set labels (subset ov ytr).
yva = bigarray_mat(fullfile(wpath,'split/yva'));
yva.Clear();
yva.Init(1000);

codes = 2*eye(1000)-1;

%% Actually generate the dataset copying data in the bigarrays.

for i = 1:numel(list)
		load(fullfile(src,list{i}));
		order = randperm(numel(num));
		stop = floor(0.2*numel(num));
		yTemp = repmat(codes(i,:),numel(num),1);

		xte = x(:,1:stop);
		xtr = x(:,stop+1:end);
		

		Yte = yTemp(1:stop,:)';
		Ytr = yTemp(stop+1:end,:)';

		% take 20% of points for each class and throw them in the validation set.

		stop = floor(0.2*size(xtr,2));
		xva = xtr(:,1:stop);
		Yva = Ytr(:,1:stop);

		fprintf('Processing file %s: %.2f%%\r',list{i}, i*100/numel(list));
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
