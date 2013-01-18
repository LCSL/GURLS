function [splits] = split_ho(X, y, opt)
% [splits] = split_ho(X, y, opt)
% Splits data into train and validation set

% INPUTS:
% -X: not used
% -y: labels matrix
% -OPT: structure of options, with the following field with default values
%       set through the defopt function:
%       -nholdouts
%       -hoproportion
% 
% OUTPUT: struct array of length OPT.NHOLDOUTS with each element gaving the following fields:
%   -tr: indices of samples to be used for training
%   -va: indices of samples to be used for testing

nSplits = opt.nholdouts;
fraction = 1-opt.hoproportion;

[n,T] = size(y);
[dummy, y] = max(y,[],2);

for t = 1:T,
	classes{t}.idx = find(y == t);
end

splits = {};
for state = 1:nSplits,
	rand(state);

%% Shuffle each class
	splits{state} = struct;
	splits{state}.tr = [];
	splits{state}.va = [];
	
	for t = 1:T,
		nSamples = numel(classes{t}.idx);
		shuffle = randperm(nSamples);
		classes{t}.idx = classes{t}.idx(shuffle);

		begidx = 1;
		endidx = floor(fraction*nSamples);

		splits{state}.tr = unique([splits{state}.tr; classes{t}.idx(begidx:endidx)]);
		splits{state}.va = unique([splits{state}.va; classes{t}.idx(endidx+1:end)]);
	end
end
