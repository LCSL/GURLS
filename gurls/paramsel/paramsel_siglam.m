function vout = paramsel_siglam(X, y, opt)
% paramsel_siglam(X, y, OPT)
% Performs parameter selection when the dual formulation of RLS is used.
% The leave-one-out approach is used.
% It selects both the regularization parameter lambda and the kernel parameter sigma.
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

[~, T]  = size(y);

if ~isprop(opt,'kernel')
    opt.newprop('kernel', struct());
	opt.kernel.type = 'rbf';
end

opt.kernel.init = 1;
opt.kernel = kernel_rbf(X,y,opt);
nsigma = numel(opt.kernel.kerrange);

PERF = zeros(opt.nsigma,opt.nlambda,T);

for i = 1:nsigma
	opt.paramsel.sigmanum = i;
	opt.kernel = kernel_rbf(X,y,opt);
	paramsel = paramsel_loocvdual(X,y,opt);
	PERF(i,:,:) = paramsel.perf;
	guesses(i,:) = paramsel.guesses;
end

% The lambda axis is redefined each time but
% it is the same for all classes as it depends
% only on K so we can still sum and minimize.
%
% We have to be a bit careful when minimizing.
%
% TODO: select a lambda for each class fixing sigma.

M = sum(PERF,3); % sum over classes

[dummy, i] = max(M(:));
[m, n] = ind2sub(size(M),i);

% opt sigma
vout.sigma = opt.kernel.kerrange(m);
% opt lambda
vout.lambdas = guesses(m,n)*ones(1,T);
