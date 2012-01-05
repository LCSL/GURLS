function = [cols] = BigColumns(bigArray, idx)

		cols = zeros(bigArray.NumItems(),numel(idx));
		mBegin = 0;
		for i = 1:bigArray.NumBlocks()
			M = bigArray.ReadBlock(i);
			M = M(idx,:);
			cols(mBegin+1:mBegin+size(M,1),:) = M';
			mBegin = mBegin + size
		end	
