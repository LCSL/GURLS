function [scores] = pred_dual_mkl(X, y, opt)
% pred_primal(X,y,opt)
% computes the predictions of the classifier stored in opt.rls on
% the samples passed in the X matrix.
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%   -X: input data matrix, not used...
%   -y: labels matrix. only used to measure n_pred
%   fields with default values set through the defopt function:
%		- kernel.type
%   fields that need to be set through previous gurls tasks:
%		- rls (set by the rls_* routines)
%       - predkernel.K (set by the predkernel_* routines), only if
%         kernel.type is not 'linear'
%
% OUTPUT:
% -scores: npred x 1 matrix of predicted labels

if strfind(opt.kernel.type, 'mkl')
    if isprop(opt, 'predkernel')
        npred = length(y);
        M = size(opt.rls.C_mkl, 2);
        scores = zeros(npred, 1);
        for m = 1:M
            scores = scores + ...
                opt.predkernel.K_mkl(:,:,m) * opt.rls.C_mkl(:, m);
        end
    else
        error('Please provide a predkernel');
    end
else
    error('opt.kernel.type = ''%s'', expect ''mkl_*''', ...
        opt.kernel.type);
end
