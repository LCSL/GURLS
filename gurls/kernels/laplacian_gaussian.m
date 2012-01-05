function [laplacian] = laplacian_gaussian(X,y,opt)


if isfield(opt,'kernel')
	if isfield(opt.kernel,'distance')
		opt.laplacian.distance = opt.kernel.distance;
	else
		opt.laplacian.distance = distance(X',X');
	end
else	
	opt.laplacian.distance = distance(X',X');
end	

%D = -(opt.kernel.distance.^2);
W = exp((-opt.laplacian.distance)/(opt.paramsel.sigma^2));
D = diag(sum(W));
L = D-W;
laplacian.L = L;

laplacian.type = 'gaussian';
