function [ K ] = preproc_kernel_dispatch( kernel_opts, X1, X2 )

    % There doesn't seem to be any good way of fitting this in the GURLS
    % framework.  It looks like whoever implemented predkernel_traintest.m
    % had to re-implement all of the kernels, so I guess I'll do it here.
    size(X1)
    size(X2)
    switch kernel_opts.kernel
        case 'linear'
            K = X1*X2';
        case 'rbf'
            D = square_distance(X1',X2'); 
            K = exp(-0.5*D/kernel_opts.sigma^2);
        case 'poly'
            K = X1*X2';
            K = (K + kernel_opts.c).^kernel_opts.order;
        otherwise
            error('Only linear, RBF, and polynomial kernels are supported for preprocessing.');
    end
end

