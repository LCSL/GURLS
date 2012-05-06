function [X] = norm_testzscore(X, y, opt)
% norm_testzscore(X,Y,OPT)
% Normalizes test data using the same mean and std deviation computed for the training set, 
% previously computed and stored in the file with name root specified in opt.name by the norm_zscore function.		
%
% INPUTS:
% -X: input data matrix
% -y: not used
% -OPT: structure of options with the following fields (and subfields):
%		- opt.name
% 
% OUTPUT:
% -X: normalized input data matrix


load([opt.name '_norm_zscore.mat' ]);

[n,d] = size(X);
X = X - repmat(meanX, n, 1);
X = X./repmat(stdX, n, 1);
