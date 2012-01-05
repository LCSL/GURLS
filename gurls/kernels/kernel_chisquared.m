function [kernel] = kernel_chisquared(X,y,opt)
		for i = 1:size(X,1)
			for j = 1:i
				kernel.K(i,j) = sum(...
								( (X(i,:) - X(j,:)).^2 ) ./ ...
								( 0.5*(X(i,:) + X(j,:)) + eps));
				kernel.K(j,i) = kernel.K(i,j);
			end
		end	
		kernel.type = 'chisquared';
