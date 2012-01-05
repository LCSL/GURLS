function [splits] = split_monline(X, y, opt)
% Splits the dataset into opt.nSplits number of splits
% such that each split has equal number of examples per class.

if isfield (opt,'split')
	splits = opt.split;	% lets not overwrite existing performance measures.
			   	      	% unless they have the same name
	if isfield (opt.split,'tr')
		y = y(opt.split.tr,:);
	end	
end

opt.nSplits = opt.nSplits + 1;

T = size(y,2);

[dummy, y] = max(y,[],2);

for t = 1:T,
	classes{t}.idx = find(y == t);
	classes{t}.splitidx = floor(linspace(0,numel(classes{t}.idx),opt.nSplits));
end
%% Shuffle each class
for t = 1:T,
	nSamples = numel(classes{t}.idx);
	shuffle = randperm(nSamples);
	classes{t}.idx = classes{t}.idx(shuffle);
end

for j = 1:opt.nSplits-1,
	bag_splits{j}.idx = [];
% 	for t = 1:10,
    for t = 1:T,
		begidx = classes{t}.splitidx(j) + 1;
		endidx = classes{t}.splitidx(j + 1);
		bag_splits{j}.idx = [bag_splits{j}.idx; classes{t}.idx(begidx:endidx)]; 
	end
	splits.bag = bag_splits;
end

