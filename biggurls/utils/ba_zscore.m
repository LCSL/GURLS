function [m,s] = ba_zscore(B, m, s)

		% If m and s are not sent, compute them.
		if (nargin < 2)
			[m,s] = ba_meanstd(B);
		end

		B.BlockOp(@whatever,1)
		function [X] = whatever(X)
			for i=1:size(X,1)
				X(i,:) = X(i,:) - m(i);
				X(i,:) = X(i,:)/(s(i)+eps);
			end	
		end
end

