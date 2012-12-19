function [bigMat] = bigLoad(bigArray)
		mBegin = 0;
		for i = 1:bigArray.NumBlocks()
			M = bigArray.ReadBlock(i);
			M = M';
			bigMat(mBegin+1:mBegin + size(M,1), :) = M;
			mBegin = mBegin + size(M,1);
		end
end		

