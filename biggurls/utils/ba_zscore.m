function [m,s] = ba_zscore(B, m, s)

		% If m and s are not sent, compute them.
		if (nargin < 2)
			[m,s] = ba_std(B);
		end

		B.BlockOp(@whatever,1)
		function [X] = whatever(X)
			for i=1:size(X,1)
				X(i,:) = X(i,:) - m(i);
				X(i,:) = X(i,:)/(s(i)+eps);
			end	
		end
end

function [m, s] = ba_std(B, prep_func)
	
	ex = 0;
	exsq = 0;
	n = 0;
	

	if nargin < 2
		prep_func = @(X) X;
	end

	m = 0;

	B.BlockOp(@block_sum);

	function block_sum(X)

		X = prep_func(X);
		n = n + size(X,2);

		ex = ex + sum(X,2);
		exsq = exsq + sum(X.^2,2);
	end

	if n > 0
		m = ex/n;
		s = sqrt(exsq - ex.^2)/n;
	else
		error('No elements match your prep_func criterion');
	end

end
