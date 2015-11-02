function vout = paramsel_siglamho_nystrom(X, y, opt)
% paramsel_siglamho_nystrom(X, y, OPT)
% Performs parameter selection when Nystrom KRLS is used.
% The hold-out approach is used.
% It selects both the regularization parameter lambda and the kernel parameter sigma.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- kernel.Knm (set by the kernel_*_nystrom routines)
%		- kernel.Kmm (set by the kernel_*_nystrom routines)
%   fields with default values set through the defopt function:
%       - nsigma
%		- kernel.type
%		- nlambda
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

[~,T]  = size(y);

if ~isprop(opt,'kernel')
    opt.newprop('kernel', struct());
	opt.kernel.type = 'rbf_nystrom';
end

opt.kernel.init = 1;
opt.kernel = kernel_rbf_nystrom(X,y,opt);
nsigma = numel(opt.kernel.kerrange);

if ~exist('vout','var') || ~isfield(vout, 'regrange')
    tot = opt.nlambda;
else
    tot = numel(vout.regrange);
end

PERF = zeros(opt.nsigma,tot,T);

for i = 1:nsigma
	opt.paramsel.sigmanum = i;
	opt.paramsel.sigma = opt.kernel.kerrange(i);
	paramsel = paramsel_hodual_nystrom(X,y,opt);
	nh = numel(paramsel.perf);
    nl = numel(paramsel.guesses{1});
	PERF(i,:,:) = reshape(median(reshape(cell2mat(paramsel.perf')',nl*T,nh),2),T,nl)';
	guesses(i,:) = median(cell2mat(paramsel.guesses'),1);
end
% The lambda axis is redefined each time but
% it is the same for all classes as it depends
% only on K so we can still sum and minimize.
%
% We have to be a bit careful when minimizing.
%
% TODO: select a lambda for each class fixing sigma.

M = sum(PERF,3); % sum over classes

[~,i] = max(M(:));
[m,n] = ind2sub(size(M),i);

% opt sigma
vout.sigma = opt.kernel.kerrange(m);
% opt lambda
vout.lambdas = guesses(m,n)*ones(1,T);
