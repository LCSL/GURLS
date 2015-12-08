function D = dist_euclid(X1,X2,A)


if ~exist('A','var')

    X1X1 = (sum(X1.^2,2));
    if ~exist('X2','var') || numel(X2)==0
	D = ones(length(X1X1),1)*X1X1' + X1X1*ones(1,length(X1X1)) - 2*(X1*X1');
    else
	X2X2 = (sum(X2.^2,2));
	D = ones(length(X1X1),1)*X2X2' + X1X1*ones(1,length(X2X2)) - 2*(X1*X2');
    end
else
    X1X1 = (sum((X1*A).^2,2));
    
    if ~exist('X2','var') || numel(X2)==0
	X2X2 = X1X1;
    else
	X2X2 = (sum((X2*A).^2,2));
    end
    
    D = ones(length(X1X1),1)*X2X2' + X1X1*ones(1,length(X2X2)) - 2*(X1*A*X2');

end

