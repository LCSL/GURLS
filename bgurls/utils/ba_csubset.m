function [bXsub] = ba_csubset(bX, bXsub, ids)
% bX : input bigarray
% bXsub : output bigarray (subset of bX, taken columnwise)
% ids   : indices on bX to be written to bXsub.
	
	nva = numel(ids);
	
	b = bX.BlockSize();
	
	bXsub.Clear();
	bXsub.Init(b);
	
	for i = 1:b:nva,
		selids = ids(i:min(i+b-1, nva));
		bXsub.Append(bX(:, selids));
	end

	bXsub.Flush();
end
