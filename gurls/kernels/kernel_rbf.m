function [kernel] = kernel_rbf(X, y, opt)

% 	kernel_rbf(opt)
%	Computes the kernel matrix for a Gaussian kernel.
%	INPUTS:
%		-OPT: struct with the following options:
%			- paramsel: struct containing the following fields (computed by paramsel_*).
%			- sigma: width of the gaussian kernel.
%           -X: input data matrix
%
%	OUTPUT: struct with the following fields:
%		-type: 'rbf'
%		-K: kernel matrix

if ~isprop(opt, 'kernel')
    opt.newprop('kernel', struct());
end
kernel = opt.kernel;

if ~isfield(kernel, 'distance')
    kernel.distance = square_distance(X', X');
end

if ~isfield(kernel, 'kerrange') && isprop(opt, 'nsigma')
    % compute kernel param range only if nsigma is set in opt
    
    if ~isprop(opt, 'sigmamin')
        n = size(kernel.distance, 1);
        D = sort(kernel.distance(tril(true(n), -1)));
        firstPercentile = round(0.01*numel(D) + 0.5);
        opt.newprop('sigmamin', sqrt(D(firstPercentile)));
        clear D n;
    end
    if ~isprop(opt, 'sigmamax')
        opt.newprop('sigmamax', sqrt(max(max(kernel.distance))));
    end
    if opt.sigmamin <= 0
        opt.sigmamin = eps;
    end
    if opt.sigmamax <= 0
        opt.sigmamax = eps;
    end

    q = (opt.sigmamax/opt.sigmamin)^(1/(opt.nsigma - 1));
    kernel.kerrange = opt.sigmamin*(q.^(opt.nsigma:-1:0));
end

if ~isfield(kernel, 'init') || ~kernel.init
    
    if ~isfield(opt.paramsel, 'sigma')
        sigma = kernel.kerrange(opt.paramsel.sigmanum);
    else
        sigma = opt.paramsel.sigma;
    end
    
    D = -(kernel.distance);
    K = exp(D/(sigma^2));
    
    kernel.type = 'rbf';
    kernel.K = K;
else
    kernel.init = 0;
end
