function [cfr] = rls_insta (X, y, opt)

% rls_insta(X,y,opt)
% computes a classifier for elastic nerwork using ISTA.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
%   fields that is optional
%       - paramsel.insta_alpha (paramters for balancing l1-norm and l-2 norm. 1 for LASSO and 0 for ridge)
%       - paramsel.niter (maximun number for iteration. Set to either negative number or inf for using threshold rule only)
%       - paramsel.relthre (relative convergence threshold for iteration to stop)
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
if isfield(opt.paramsel, 'insta_alpha')
    insta_alpha=opt.paramsel.insta_alpha;
    if insta_alpha <= 0 || insta_alpha > 1
        error('Invalid alpha');
    end
else
    if opt.verbose
            warning('alpha not defined. Use default value alpha=1 for LASSO');
            insta_alpha=1;
    end
end

% load in number of iterations or relative tolerance
Niter=-1;
relthre=1e-4;
if isfield(opt.paramsel, 'niter')
    Niter=opt.paramsel.niter;
end
if isfield(opt.paramsel, 'relthre')
    relthre=opt.paramsel.relthre;
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


% redo OLR based on non-sparsy components
% w = rls_insta_driver( XtX, Xty, n, lambda,insta_alpha,Niter,relthre,opt);
% cfr.IndexFlag = ~~(w);
% Xs=X(:,~~w);
% cfr.W=zeros(size(w));
% cfr.W(~~w) = rls_primal_driver(Xs'*Xs,Xs'*y,n,0);
cfr.W = rls_insta_driver( XtX, Xty, n, lambda,insta_alpha,Niter,relthre,opt);
cfr.C = [];
cfr.X = [];

