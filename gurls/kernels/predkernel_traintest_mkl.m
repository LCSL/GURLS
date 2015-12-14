function [fk] = predkernel_traintest_mkl(X, y, opt)
% predkernel_traintest(OPT)
% Computes the kernel matrix between the training points and the
% test points. It can be used by pred_dual.
%
% INPUTS:
% -OPT: structure with the following fields:
%   fields that need to be set through previous gurls tasks:
%       - X: input testing data matrix
%		- rls.X: input training data matrix
%   fields (mkl-related) that need to be set by hand:
%       - mkl.type: kernel type and specification
%                   (see kernel_mkl for more)
%   For more information on standard OPT fields
%   see also defopt
%
% OUTPUT: struct (opt.predkernel) with at least
%         the field K_mkl containing the ntest x ntr x M kernel matrix

if ~isprop(opt,'predkernel')
    opt.newprop('predkernel', struct());
    opt.predkernel.type = 'mkl';
end

X_tr = opt.rls.X;
X_va = X;
n_tr = size(opt.rls.X, 1);
n_va = size(X, 1);
M = size(opt.kernel.K_mkl, 3);

mkl_type = opt.mkl.type;
K_list = zeros(n_va, n_tr, M);
K_id = 0;

for idx = 1:length(mkl_type)
    spec = mkl_type{idx};
    kern_fun = str2func(spec{1});
    par_list = spec{2};
    % evaluate kernels: rbf case
    if strcmp(spec{1}, 'kernel_rbf')
        rbf_dist = square_distance(X_va', X_tr');
        for sigma_idx = 1:length(par_list)
            K_id = K_id + 1;
            K_temp = exp(-rbf_dist/par_list(sigma_idx));
            K_list(:, :, K_id) = K_temp;
        end
    elseif strcmp(spec{1}, 'kernel_linear')
        K_id = K_id + 1;
        K_temp = X_va * X_tr';
        K_list(:, :, K_id) = K_temp;
    else
        fprintf('''%s'' currently not supported, skip\n', spec{1})
    end
end

fk.K_mkl = K_list;

end


