function vout = paramsel_siglamloogpregr(X,y, opt)
% paramsel_siglamloogpregr(X,y, OPT)
% Performs parameter selection for gaussian process regression.
% The leave-one-out approach is used.
% It selects both the noise level lambda and the kernel parameter sigma.
%
% INPUTS:
% -OPT: struct of options with the following fields:
%   fields with default values set through the defopt function:
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

%savevars = {'LOOSQE','M','sigmas','guesses'};

if isprop(opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
else
    opt.newprop('paramsel', struct());
end

[~,T]  = size(y);

if ~isprop(opt,'kernel')
    opt.newprop('kernel', struct());
	opt.kernel.type = 'rbf';
end

kerfun = str2func(['kernel_' opt.kernel.type]);
opt.kernel.init = 1;
opt.kernel = kerfun(X,y,opt);
nsigma = numel(opt.kernel.kerrange);


for i = 1:nsigma
	opt.paramsel.sigmanum = i;
	opt.kernel = kerfun(X,[],opt);
	paramsel = paramsel_loogpregr(X,y,opt);
	perf(i,:,:) = paramsel.perf;
	guesses(i,:) = paramsel.guesses;
end
% The lambda axis is redefined each time but
% it is the same for all classes as it depends
% only on K so we can still sum and minimize.
%
% We have to be a bit careful when minimizing.
%
% TODO: select a lambda for each class fixing sigma.

M = sum(perf,3); % sum over classes

[dummy,i] = max(M(:));
[m,n] = ind2sub(size(M),i);

vout.sigma = opt.kernel.kerrange(m);
vout.lambdas = guesses(m,n)*ones(1,T);

if numel(savevars) > 0
	[ST,I] = dbstack();
	save(ST(1).name,savevars{:});
end	
