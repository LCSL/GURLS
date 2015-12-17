function vout = paramsel_siglam_semi(X, y, opt)
% paramsel_semisup(X, y, OPT)
% Performs parameter selection for semi-supervised learning with the dual 
% formulation of RLS and manifold regularization.
% The leave-one-out approach is used.
% It selects both the regularization parameter lambda_m and the kernel
% parameter sigma_m for the manifold regularization, and the regularization 
% paramater lambda and the kernel parameter sigma for the ambient 
% regularization term.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- kernel.K (set by the kernel_* routines)
%   fields with default values set through the defopt function:
%		- kernel.type
%		- nlambda
%		- nsigma
%       - hoperf
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: structure with the following fields:
% -lambdas: value of the regularization parameter lambda
%           minimizing the validation error, replicated in a TX1 array 
%           where T is the number of classes
% -sigma: value of the kernel parameter minimizing the validation error
% -lambdas_m: value of the regularization parameter lambda
%             minimizing the validation error, replicated in a TX1 array 
%             where T is the number of classes
% -sigma_m: value of the kernel parameter minimizing the validation error

%savevars = {'LOOSQE','M','sigmas','guesses'};

if isprop(opt,'paramsel')
    vout = opt.paramsel; % lets not overwrite existing parameters.
    % unless they have the same name
else
    opt.newprop('paramsel', struct());
end

% case: sigma was prev. set & paramsel_siglam is called to re-compute it
if isfield(opt.paramsel, 'sigma')
    opt.paramsel = rmfield(opt.paramsel, 'sigma');
end

[n,T] = size(y);

if ~isprop(opt,'kernel')
    opt.newprop('kernel', struct());
	opt.kernel.type = 'rbf';
end

opt.kernel.init = 1;
opt.kernel = kernel_rbf(X,y,opt);

% find range of labelled y
labelled = (~isnan(y(:,1)));
labelledrange = find(labelled);
J = diag(+labelled);
m = sum(+labelled);
Jy = zeros(size(y));
Jy(labelled,:) = y(labelled,:);
D = -(opt.kernel.distance);
                
nsigma = numel(opt.kernel.kerrange);
nlambda = opt.nlambda;
PERF = zeros(nsigma,nsigma,nlambda,nlambda,T);
lambda_guesses = zeros(nsigma,nlambda);


for i = 1:nsigma
    % find the guesses_m
    opt.paramsel.sigmanum = i;
    opt.kernel = kernel_rbf(X,y,opt);
    [~,L] = eig(opt.kernel.K);
    L_m = double(diag(L));
    guesses_m = paramsel_lambdaguesses(L_m, n, n, opt);
    sigma_m = opt.kernel.kerrange(i);
    lambda_guesses(i,:) = guesses_m;

    for j = 1:nsigma 
        % find the guesses
        opt.paramsel.sigmanum = j;
        opt.kernel = kernel_rbf(X,y,opt);
        [~,L] = eig(opt.kernel.K);
        L = double(diag(L));
        guesses = paramsel_lambdaguesses(L, n, n, opt);

        for k = 1:nlambda;
            for l = 1:nlambda;
                lambda_m = guesses_m(k);
                lambda = guesses(l);
                % pred using manfold_learning and find performance
                for lo = labelledrange'
                    r = [1:lo-1 lo+1:n];
                    
                    W = exp(D(r,r)./(sigma_m^2)).*(opt.kernel.distance(r,r) < 4.*sigma_m^2);
                    G = diag(sum(W,2))-W;
                    for t = 1:T
                        C =(J(r,r)*opt.kernel.K(r,r) + (m*lambda)*eye(n-1) + (m*lambda_m)*G*opt.kernel.K(r,r)) \ Jy(r,t);
                        y_pred = opt.kernel.K(lo,r)*C;
                        PERF(i,j,k,l,t) = PERF(i,j,k,l,t) - (y_pred - y(lo,t))^2;
                    end
                end
            end
        end
    end
end

M = sum(PERF,5); % sum over classes

[~, i] = max(M(:));
[mi, mj, mk, ml] = ind2sub(size(M),i);

vout.sigma_m = opt.kernel.kerrange(mi);
vout.sigma = opt.kernel.kerrange(mj);
vout.lambdas_m = lambda_guesses(mi,mk)*ones(1,T);
vout.lambdas = lambda_guesses(mj,ml)*ones(1,T);



