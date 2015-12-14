function w=rls_insta_driver( XtX, Xty, n, lambda,inst_alpha,Niter,relthre,opt)
% rls_insta_driver( XtX, Xty, n, lambda,alpha,Niter,reltol,opt)
% Utility function used by rls_insta
% 
% INPUTS:
% -XtX: symmetric dxd square matrix
% -Xty: dxT matrix
% -n: number of training samples
% -lambda: regularization parameter
% -alpha: l1 norm l2 norm weight parameter
% -Niter: maximum number of iteration
% -retol: relative tolernce for convergence
% -opt: use opt.verbose only
%
% OUTPUTS:
% -W: matrix of coefficient vector for linear RLS classifier

    gamma = n/(2*eigs(XtX,1)+lambda*(1-inst_alpha)); %Step size
    
    % make sure Niter, retol are proper
    if ~isinteger(Niter)
        if opt.verbose
        warning('Interation number must be integer. Rounding opt.paramsel.niter');
        end
        Niter = ceil(Niter);
    end
    if Niter<=0 && relthre<0
        if opt.verbose
            warning(['Unvalid stopping rule for INSTA.',...
            'Using default relative tolerance = 1e-4']);
        end
        Niter = -1;
        relthre=1e-4;      
    end
    
    % start insta

    w2=XtX\Xty;
    w2(abs(w2)<lambda*inst_alpha*gamma)=0;
    w2=w2-sign(w2)*lambda*inst_alpha*gamma;
    w1=0*w2;
    %dw = mean(abs(w1-w2));
    k=0;
    while mean(abs(w1-w2))>relthre*max([mean(abs(w1)),mean(abs(w2))])
        w1=w2;
        w2=(1-lambda*gamma*(1-inst_alpha))*w2+2*gamma*(Xty-XtX*w2)/n;
        w2(abs(w2)<lambda*inst_alpha*gamma)=0;
        plot(w2)
        w2=w2-sign(w2)*lambda*inst_alpha*gamma;
        k=k+1;

        if k==Niter
            break;
        end

    end
    w=w2;
end