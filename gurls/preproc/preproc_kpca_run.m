function [preproc,X,y] = preproc_kpca_run(X,y,opt)
    preproc = opt.preproc;
    
    if strcmp(preproc.kernel.kernel,'linear')
        % special case: this makes PCA (with no 'K') faster
        X = X*preproc.V;
    else
        K = preproc_kernel_dispatch(preproc.kernel,preproc.X,X);
        X = K*preproc.V;
    end
end