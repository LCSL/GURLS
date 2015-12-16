function [kernel] = kernel_mkl(X, y, opt)

% 	kernel_rbf(opt)
%	Computes M n x n kernel matrices given input
%	INPUTS:
%		-OPT: struct with the following options:
%			- mkl.type: cell array with element {'type', params}
%             support: type = rbf: params = 1 x M_rbf array of sigmas
%                      type = linear: params = 0
%       -X: input data matrix
%
%	OUTPUT: struct with the following fields:
%		-type: 'mkl'
%		-K_mkl: n x n x M array of M n x n kernel matrices
%		-eig_mkl: M x 1 array of max eigenvalues of M kernel matrices

% initialize output container
if ~isprop(opt,'kernel')
    opt.newprop('kernel', struct());
end
mkl_type = opt.mkl.type;
kernel = opt.kernel;

% compute kernel specified by mkl_type
if any(cellfun('length', mkl_type)~=2)
    error('parameter not supplied for some kernels')
end

n = size(X, 1);
M = sum(cellfun(@(x) length(x{2}), mkl_type));
K_list = zeros(n, n, M);
eig_list = zeros(M, 1);
K_id = 0;

for idx = 1:length(mkl_type)
    spec = mkl_type{idx};
    kern_fun = str2func(spec{1});
    par_list = spec{2};
    % evaluate kernels: rbf case
    if strcmp(spec{1}, 'kernel_rbf')
        for sigma_idx = 1:length(par_list)
            K_id = K_id + 1;
            opt = set_sigma(opt, par_list(sigma_idx));
            K_temp = kern_fun(X, y, opt);
            K_list(:, :, K_id) = K_temp.K;
            eig_list(K_id) = eigs(K_temp.K, 1);
        end
    elseif strcmp(spec{1}, 'kernel_linear')
        K_id = K_id + 1;
        K_temp = kern_fun(X, y, opt);
        K_list(:, :, K_id) = K_temp.K;
        eig_list(K_id) = eigs(K_temp.K, 1);
    else
        fprintf('''%s'' currently not supported, skip\n', spec{1})
    end
end

kernel.type = 'mkl';
kernel.K_mkl = K_list;
kernel.eig_mkl = eig_list;

if isprop(opt, 'paramsel')
    if isfield(opt.paramsel, 'manual_sigma')
        % clear opt.paramsel if it is defined by @set_sigma
        opt.paramsel = struct();
    end
end
end

function [opt] = set_sigma(opt, sigma)
% insert sigma to the slot 'opt.paramsel.sigma'
if ~isprop(opt,'paramsel')
    opt.newprop('paramsel', struct());
    % signature
    opt.paramsel.manual_sigma = true;
end
opt.paramsel.sigma = sigma;
end

