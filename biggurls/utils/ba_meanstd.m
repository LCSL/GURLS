function [m,s] = ba_meanstd(B, prep_func)

	ex = zeros(B.Sizes(1),1);
	exsq = zeros(B.Sizes(1),1);
	n = zeros(B.Sizes(1),1);
	

	if nargin < 2
		prep_func = @(X) X;
	end

	m = 0;

	B.BlockOp(@block_sum);

	function block_sum(X)
		X = prep_func(X);
		sumX = zeros(size(X,1),1);
		numX = zeros(size(X,1),1);
		sumXSQ = zeros(size(X,1),1);
		for i = 1:size(X,1)
			values = X(i,:);
			sumX(i) = sum(values);
			sumXSQ(i) = sum(values.^2);
			numX(i) = numel(values);
		end	
		
		ex = ex + sumX;
		exsq = exsq + sumXSQ;
		n = n + numX;
	end
	if n > 0
		m = ex./n;
		s = sqrt(exsq./n - (ex./n).^2);
	else
		error('No elements match your prep_func criterion');
	end

end
