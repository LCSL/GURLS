function [splits] = splits_trva(X, y, opt, nSplits, fraction)
% [splits] = splits_incbatch(X, y, opt)
% Expects binary coded 'y', just like all other functions

if nargin < 4, nSplits = 10; end
if nargin < 5, fraction = 0.8; end

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
