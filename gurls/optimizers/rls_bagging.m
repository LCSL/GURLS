function [cfr] = rls_bagging(X, y, opt)
% Split X into opt.nSplits sets
nSplits = opt.nSplits;
disp(sprintf('- number of splits = %d',nSplits));

using_bigarray = 0;
%if strcmp(class(X),'bigarray') || strcmp(class(X),'mmapbigarray')
if isa(X, 'bigarray')
	using_bigarray = 1;
end

if using_bigarray
	n = X.NumItems();
else
	[n, d] = size(X);
end
T = size(y,2);

cfr = struct;
try
	opt.expname = 'bagging_dummy';
	opt.savefile = ['bagging_dummy.mat'];
	opt.seq = opt.bagseq;
	opt.process = opt.bagproc;
catch
	error('Specify opt.bagseq and opt.bagproc for training individual models\n');
end

bag_splits = make_splits(y, nSplits+1,'macrobatch');

for splitno = 1:nSplits,
	fprintf('\n\t Split [%d]',splitno);

	opt.singlelambda = @median;
	sel_ids = bag_splits{splitno}.idx;

	fprintf('%s\n',opt.expname);

	if using_bigarray
		sel_ids = sort(sel_ids);
		M = X(sel_ids,:);
		gurls(M, y(sel_ids,:), opt,1);
	else
		gurls(X(sel_ids,:), y(sel_ids,:), opt,1);
	end

	t = load([opt.savefile]);
	fprintf('Loaded %s\n',t.opt.savefile);

	cfr.W{splitno} = t.opt.rls.W;
	cfr.C{splitno} = t.opt.rls.C;
	cfr.X{splitno} = t.opt.rls.X;
	fprintf('\t done..');
end
