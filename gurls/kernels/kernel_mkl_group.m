function [kernel] = kernel_mkl_group(X, y, opt)

% 	kernel_rbf(opt)
%	Computes M n x n kernel matrices to use
%       different kernels for different features
%
%	INPUTS:
%		-OPT: struct with the following options:
%			- mkl.group.type:
%             > cell array with element {'type', params}
%             support: type = rbf: params = 1 x M_rbf array of sigmas
%                      type = linear: params = 0
%           - mkl.group.idx:
%             > cell array indicating the index of feature group
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

mkl_type = opt.mkl.group.type;
mkl_group = opt.mkl.group.idx;
kernel = opt.kernel;

% compute kernel specified by mkl_type
if any(cellfun(@(x) length(x{1}), mkl_type)~=2)
    error('parameter not supplied for some kernels')
end

M_group = zeros(1, length(mkl_type));
group_id = [];
for idx = 1:length(mkl_type)
    M_group(idx) = sum(cellfun(@(x) length(x{2}), mkl_type{idx}));
    group_id = [group_id ones(1, M_group(idx))*idx];
end

M = sum(M_group);
n = size(X, 1);
K_list = zeros(n, n, M);
eig_list = zeros(M, 1);

K_id = 0;

for group_idx = 1:length(mkl_group)
    X_group = X(:, mkl_group{group_idx});
    mkl_type_group = mkl_type{group_idx};
    for idx = 1:length(mkl_type_group)
        spec = mkl_type_group{idx};
        kern_fun = str2func(spec{1});
        par_list = spec{2};
        % evaluate kernels: rbf case
        if strcmp(spec{1}, 'kernel_rbf')
            for sigma_idx = 1:length(par_list)
                K_id = K_id + 1;
                opt = set_sigma(opt, par_list(sigma_idx));
                K_temp = kern_fun(X_group, y, opt);
                K_list(:, :, K_id) = K_temp.K;
                eig_list(K_id) = eigs(K_temp.K, 1);
            end
        elseif strcmp(spec{1}, 'kernel_linear')
            K_id = K_id + 1;
            K_temp = kern_fun(X_group, y, opt);
            K_list(:, :, K_id) = K_temp.K;
            eig_list(K_id) = eigs(K_temp.K, 1);
        else
            fprintf('''%s'' currently not supported, skip\n', spec{1})
        end
    end
end

if opt.mkl.group.interaction
    if length(mkl_group) == 2
        num_int_kernel = prod(M_group);
        user_agree = ...
            input(sprintf(...
            '\n %d interaction kernels will be created, continue?(y/n)', ...
            num_int_kernel), 's');
        if strcmp(user_agree, 'y')
            % create interaction kernels after user consent
            % 1. container and index list for kernel pairs
            Kint_list = zeros(n, n, num_int_kernel); % container
            eigint_list = zeros(M, 1);

            Kid_list = mat2cell(1:M, 1, M_group);
            [p,q] = meshgrid(Kid_list{1}, Kid_list{2});
            Kint_pair = [p(:) q(:)]; % index list
            
            % 2. interaction kernel creation
            for Kint_id = 1:num_int_kernel
                id_pair = Kint_pair(Kint_id, :);
                K_temp = ...
                    K_list(:,:,id_pair(1)).* K_list(:,:,id_pair(2));
                Kint_list(:, :, Kint_id) = K_temp;
                eigint_list(Kint_id) = eigs(K_temp, 1);
            end
            
            % 3. merge back to orginal K_list/eig_list
            K_list = cat(3, K_list, Kint_list); % container
            eig_list = [eig_list; eigint_list];
            group_id = [group_id, 3*ones(1, num_int_kernel)];
        else 
            fprintf('no interaction kernel created\n');
        end
    else
       error('Expect 2 feature groups, %d given', length(mkl_group));
    end
end

kernel.type = 'mkl_group';
kernel.K_mkl = K_list;
kernel.eig_mkl = eig_list;
kernel.group_idx = group_id;

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