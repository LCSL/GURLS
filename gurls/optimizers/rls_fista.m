function [cfr] = rls_fista (X, y, opt)
% rls_fista(X,y,opt)
% computes a classifier for elastic nerwork using FISTA.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
%   fields with default values set through the defopt function:
%       - ISTAalpha (paramters for balancing l1-norm and l-2 norm. 1 for LASSO and 0 for ridge)
%       - niter (maximun number for iteration. Set to either negative number or inf for using threshold rule only)
%       - relthre (relative convergence threshold for iteration to stop)
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -C: matrix of coefficient vectors of rls estimator for each class
% -X: empty matrix

lambda = opt.singlelambda(opt.paramsel.lambdas);

n = size(y,1);

% load in parameters alpha
if isprop(opt,'ISTAalpha')
    ISTAalpha=opt.ISTAalpha;
    if ISTAalpha <= 0 || ISTAalpha > 1
        error('Invalid alpha');
    end
else
    if opt.verbose
            fprintf('\t...alpha not defined. Use default value alpha=1 for LASSO\n');
            ISTAalpha=1;
    end
end

% load in number of iterations or relative tolerance
Niter=inf;
relthre=1e-4;
if isprop(opt, 'ISTAniter')
    Niter=opt.ISTAniter;
end
if isprop(opt, 'ISTArelthre')
    relthre=opt.ISTArelthre;
end

% check if matrices XtX and Xty have been previously computed during
% parameter selection
if isfield(opt.paramsel,'XtX');
    XtX = opt.paramsel.XtX;
else
    XtX = X'*X; % d x d matrix.
end
if isfield(opt.paramsel,'XtX');
    Xty = opt.paramsel.Xty;
else
    Xty = X'*y; % d x T matrix.
end


%redo OLR based on non-sparsy components


w = rls_fista_driver( XtX, Xty, n, lambda,ISTAalpha,Niter,relthre,opt.verbose);

cfr.IndexFlag = ~~(w);
% Xs=X(:,~~w);
% cfr.Wr=zeros(size(w));
% cfr.Wr(~~w) = rls_primal_driver(Xs'*Xs,Xs'*y,n,0);
cfr.W = w;
cfr.C = [];
cfr.X = [];
%plot(w)

