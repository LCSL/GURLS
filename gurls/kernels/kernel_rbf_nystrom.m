function [kernel, nystrom] = kernel_rbf_nystrom(X, y, opt)

% 	kernel_rbf(opt)
%	Computes the Knm and Kmm Nystrom kernel approximation matrices for a Gaussian kernel.
%	INPUTS:
%		-OPT: struct with the following options:
%			- paramsel : struct containing the following fields (computed by paramsel_*).
%				- sigma : width of the gaussian kernel.
%           -X: input data matrix
%           -nystrom.m: Nystrom subsampling level (a.k.a. number of sampled
%           columns of K)
%           -nystrom.sampledIdx: list of m sampled indexes of training
%           points to be used for Knm, Kmm computation after model
%           selection
%	
%	OUTPUT: struct with the following fields:
%		-type: 'rbf_nystrom'
%		-kernel.Knm: m sampled columns of the full kernel matrix
%		-kernel.Kmm: m sampled columns and rows of the full kernel matrix

if ~isprop(opt,'kernel')
    opt.newprop('kernel', struct());
end
kernel = opt.kernel;

ntr = size(X,1);

% Compute sigma range

if ~isfield(kernel, 'kerrange')
    
    sq_dist_mat = square_distance(X(1:min([ntr,2000]),:)',X(1:min([ntr,2000]) , :)');
    
    if ~isprop(opt,'sigmamin')
        D = sort(sq_dist_mat(tril(true(ntr),-1)));
        firstPercentile = round(0.01*numel(D)+0.5);
        opt.newprop('sigmamin', sqrt(D(firstPercentile)));
        clear D;
    end
    if ~isprop(opt,'sigmamax')
        opt.newprop('sigmamax', sqrt(max(max(sq_dist_mat))));
    end
    if opt.sigmamin <= 0
        opt.sigmamin = eps;
    end
    if opt.sigmamax <= 0
        opt.sigmamax = eps;
    end	
    q = (opt.sigmamax/opt.sigmamin)^(1/(opt.nsigma-1));
    kernel.kerrange = opt.sigmamin*(q.^(opt.nsigma:-1:0));

    clear sq_dist_mat;
end


if ~isfield(kernel, 'init') || ~kernel.init
    
    % Sample m columns and compute
    nystrom.sampledIdx = 1:opt.nystrom.m;
    
    % Compute Knm, Kmm
    if ~isfield(kernel,'distance')
        kernel.distance = square_distance(X(:,:)',X(nystrom.sampledIdx,:)');
        
    elseif isfield(kernel,'distance') && (size(X,1) ~= size(kernel.distance,1))
        
        % Kernel recomputation after cross validation has been completed
        nystrom.sampledIdx = opt.nystrom.sampledIdx;
        kernel.distance = square_distance(X(:,:)',X(nystrom.sampledIdx,:)');
    end

    if ~isfield(opt.paramsel, 'sigma')
        sigma = kernel.kerrange(opt.paramsel.sigmanum);
    else
        sigma = opt.paramsel.sigma;
    end
    
    D = -(kernel.distance);
    kernel.Knm = exp(D/(sigma^2));
    kernel.Kmm = kernel.Knm(nystrom.sampledIdx , :);

    kernel.type = 'rbf_nystrom';
else
    kernel.init = 0;
end
