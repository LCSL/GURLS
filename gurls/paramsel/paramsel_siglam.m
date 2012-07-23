function vout = paramsel_siglam(X,y,opt)
% paramsel_siglam(X,Y,OPT)
% Performs parameter selection when the dual formulation of RLS is used.
% The leave-one-out approach is used.
% It selects both the regularization parameter lambda and the kernel parameter sigma.
%
% INPUTS:
% -X: input data matrix
% -Y: labels matrix
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- kernel.K (set by the kernel_* routines)
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


if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
end
savevars = [];
[n,T]  = size(y);
if ~isfield(opt,'kernel')
	opt.kernel.type = 'rbf';
end
if ~isfield(opt.kernel,'distance')
	opt.kernel.distance = distance(X',X');
end	
if ~isfield(opt,'sigmamin')
	%D = sort(opt.kernel.distance);
	%opt.sigmamin = median(D(2,:));
	D = sqrt(sort(squareform(opt.kernel.distance)));
	firstPercentile = round(0.01*numel(D)+0.5);
	opt.sigmamin = D(firstPercentile);
end
if ~isfield(opt,'sigmamax')
	%D = sort(opt.kernel.distance);
	%opt.sigmamax = median(D(n,:));
	opt.sigmamax = sqrt(max(max(opt.kernel.distance)));
end
if opt.sigmamin <= 0
	opt.sigmamin = eps;
end
if opt.sigmamin <= 0
	opt.sigmamax = eps;
end	
q = (opt.sigmamax/opt.sigmamin)^(1/(opt.nsigma-1));
PERF = zeros(opt.nsigma,opt.nlambda,T);
sigmas = zeros(1,opt.nsigma);

for i = 1:opt.nsigma
	sigmas(i) = (opt.sigmamin*(q^(i-1)));
	opt.paramsel.sigma = sigmas(i);
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

[dummy,i] = max(M(:));
[m,n] = ind2sub(size(M),i);

% opt sigma
vout.sigma = opt.sigmamin*(q^(m-1));
% opt lambda
vout.lambdas = guesses(m,n)*ones(1,T);
