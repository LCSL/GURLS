function [splits] = make_splits(y, nSplits, type, state)

T = size(y,2);

[dummy, y] = max(y,[],2);

for t = 1:T,
	classes{t}.idx = find(y == t);
	classes{t}.splitidx = floor(linspace(0,numel(classes{t}.idx),nSplits));
end
%% Shuffle each class
if (nargin > 3)
	rand(state);
end
for t = 1:T,
	nSamples = numel(classes{t}.idx);
	shuffle = randperm(nSamples);
	classes{t}.idx = classes{t}.idx(shuffle);
end


switch type
	case {'incbatch'}
%% Combined splits for batch
for j = 1:nSplits-1,
	batch_splits{j}.idx = [];
	for t = 1:T
		begidx = 1; 
		endidx = classes{t}.splitidx(j + 1);
		batch_splits{j}.idx = [batch_splits{j}.idx; classes{t}.idx(begidx:endidx)]; 
	end
end
	splits = batch_splits;

	case {'macrobatch'}
%% Separate splits for bagging
for j = 1:nSplits-1,
	bag_splits{j}.idx = [];
% 	for t = 1:10,
    for t = 1:T,
		begidx = classes{t}.splitidx(j) + 1;
		endidx = classes{t}.splitidx(j + 1);
		bag_splits{j}.idx = [bag_splits{j}.idx; classes{t}.idx(begidx:endidx)]; 
	end
	splits = bag_splits;
end
end
