function [p] = MatByBig(left, right, opt)
		d = right.Sizes();
		d = d{1};
		p = zeros(size(left,1), d);
		mBegin = 0;
		for i = 1:right.NumBlocks
			M = right.ReadBlock(i);
			M = M';
			p = p + left*M;
		end	
end		


