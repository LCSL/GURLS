%  [tr_split,te_split] = create_cvsplit(Y,nFolds)
%
% Creates a Cross Validation split for the labels Y
% In: 
%  'Y' - labels
%  'nFoldsK' - how many splits to create
%
% out: 
%  'tr_split','te_split' - cell array of size [nFolds,1] with indices of
%  points in training and test splits
%
% Peter Gehler, August 27th 2008
% pgehler@tuebingen.mpg.de
function [tr_split,te_split] = create_cvsplit(Y,nFolds)


for f=1:nFolds
    tr_split{f} = [];
    te_split{f} = [];
end

classes = sort(unique(Y));

for c=1:numel(classes)
    nExamples = sum(Y==classes(c));
    
    if nExamples < nFolds
	error('less examples than folds');
    end
    
    nPerFold = floor(nExamples/nFolds);
    
    ind = find(Y==classes(c));
    ind = ind(randperm(numel(ind)));
    
    for f=1:nFolds-1
	te_ind = ind( (f-1)*nPerFold + (1:nPerFold))';
	tr_ind = setdiff(ind,te_ind);
	
	if size(te_ind,1) ==1, te_ind = te_ind'; end
	if size(tr_ind,1) ==1, tr_ind = tr_ind'; end

	tr_split{f} = [tr_split{f};tr_ind];
	te_split{f} = [te_split{f};te_ind];
    end
    
    te_ind = ind( (nFolds-1)*nPerFold+1:end)';
    tr_ind = setdiff(ind,te_ind);

    if size(te_ind,1) ==1, te_ind = te_ind'; end
    if size(tr_ind,1) ==1, tr_ind = tr_ind'; end
    
    tr_split{nFolds} = [tr_split{nFolds};tr_ind];
    te_split{nFolds} = [te_split{nFolds};te_ind];

    if (all(size(tr_split{nFolds})~=1))
	error('dimension error');
    end	
    if (all(size(te_split{nFolds})~=1))
	error('dimension error');
    end	
end

for i=1:nFolds
    if size(tr_split{i},1)==1
	tr_split{i} = tr_split{i}';
	assert(numel(tr_split{i})<numel(Y));
    end
    if size(te_split{i},1)==1
	te_split{i} = te_split{i}';
	assert(numel(te_split{i})<numel(Y));
    end
end

