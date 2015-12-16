function [preproc,X,y] = preproc_kpca_train(X,y,opt)
    [N,d] = size(X);
    preproc = opt.preproc;
   
    if strcmp(preproc.kernel.kernel,'linear')
        if preproc.center_data
            mu = mean(X,1); % average along samples
            X2 = bsxfun(@minus,X,mu);
            K = X2*X2';
        else
            K = X'*X;
        end
        preproc.V = do_eig(K);
    else
        % Kernel matrix
        K = preproc_kernel_dispatch(preproc.kernel,X,X);

        % Kernel centering
        if preproc.center_data
            ONE = ones(N,N)/N;
            K = K - ONE*K - K*ONE + ONE*K*ONE;
        end
        K = K/trace(K); % Normalize eigenvalues


        preproc.V = do_eig(K);
        preproc.X = X; % Required for computing the kernel on other data
    end
    
    function V=do_eig(K)
         % Compute top eigvals. Give hints to MATLABs solver to speed
         % things up a bit.
        eig_opts = struct;
        eig_opts.isreal = true;
        eig_opts.issym = true;
        if preproc.n_dims > size(K,1)
            preproc.n_dims = size(K,1)
        end
        [V,~] = eigs(K,preproc.n_dims,'LM',eig_opts);
    end
end