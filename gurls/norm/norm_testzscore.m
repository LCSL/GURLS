function [X] = norm_testzscore(X, y, opt)

%	norm_testzscore(X,y,opt)
%	Normalizes test data using the same mean and std deviation computed for the training set.		
%
%	NEEDS:
%		- opt.name

load([opt.name '_norm_zscore.mat' ]);

[n,d] = size(X);
X = X - repmat(meanX, n, 1);
X = X./repmat(stdX, n, 1);
