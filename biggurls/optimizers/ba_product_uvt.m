function [UVt] = ba_product_uvt(bU, bV)
% Does O = U*V;
% Currently assumes bU and bV has the
% same blocksize.
% Bigarrays store examples columnwise.

if (bU.BlockSize() ~= bV.BlockSize())
	error('big arrays with different blocksizes are not supported in this version.');
end

if (bU.NumBlocks() ~= bV.NumBlocks())
	error('big arrays with different number of blocksizes are not supported in this version.');
end


UVt = 0;
for i = 1:bU.NumBlocks()

	u = bU.ReadBlock(i);
	v = bV.ReadBlock(i);

	UVt = UVt + u*v';
end

