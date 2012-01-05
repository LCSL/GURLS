function [] = bgTrainPrepare(wpath)
	X = bigarray.Obj(fullfile(wpath,'big/trainX'));
	y = bigarray.Obj(fullfile(wpath,'big/trainY'));

	admMatMultPrepare(X,X,fullfile(wpath,'XtX.mat'));
	admMatMultPrepare(X,y,fullfile(wpath,'Xty.mat')); 

	Xva = bigarray.Obj(fullfile(wpath,'split/Xva'));
	yva = bigarray.Obj(fullfile(wpath,'split/yva'));

	admMatMultPrepare(Xva,Xva,fullfile(wpath,'XvatXva.mat'));
	admMatMultPrepare(Xva,yva,fullfile(wpath,'Xvatyva.mat'));
