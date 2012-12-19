function [p] = BigByMat(left,right);
	
		p = zeros(left.NumItems(),size(right,2));
		mBegin = 0;
		for i = 1:left.NumBlocks()
			M = left.ReadBlock(i);
			M = M';
			p(mBegin + 1:mBegin + size(M,1),:) = M*right(:,mBegin+1:mBegin+size(M,1));
			mBegin = mBegin + size(M,1);
		end	
end		

