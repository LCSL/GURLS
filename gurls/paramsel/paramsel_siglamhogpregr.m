function vout = paramsel_siglamhogpregr(X,y,opt)
% paramsel_siglamhogpregr(X,Y,OPT)
% Performs parameter selection for gaussian process regression.
% The leave-one-out approach is used.
% It selects both the noise level lambda and the kernel parameter sigma.
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
savevars = [];

if isfield (opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
end

[n,T]  = size(y);
if ~isfield(opt,'kernel')
	opt.kernel.type = 'rbf';
end
if ~isfield(opt.kernel,'distance')
    opt.kernel.distance = square_distance(X',X');
end	
if ~isfield(opt,'sigmamin')
	D = sort(opt.kernel.distance(tril(true(n),-1)));
	firstPercentile = round(0.01*numel(D)+0.5);
	opt.sigmamin = sqrt(D(firstPercentile));
end
if ~isfield(opt,'sigmamax')
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
	paramsel = paramsel_hogpregr(X,y,opt);
    nh = numel(paramsel.perf);
	PERF(i,:,:) = reshape(median(reshape(cell2mat(paramsel.perf')',opt.nlambda*T,nh),2),T,opt.nlambda)';
	guesses(i,:) = median(cell2mat(paramsel.guesses'),1);
end
% The lambda axis is redefined each time but
% it is the same for all classes as it depends
% only on K so we can still sum and minimize.
%
% We have to be a bit careful when minimizing.
%
% TODO: select a lambda for each class fixing sigma.

vout.sigmas = sigmas;
vout.lambda_guesses = guesses;

M = sum(PERF,3); % sum over classes
[dummy,i] = max(M(:));
[m,n] = ind2sub(size(M),i);
% opt sigma
vout.perf = M;
vout.sigma = opt.sigmamin*(q^(m-1));
% opt lambda
vout.lambdas = guesses(m,n)*ones(1,T);
% This is awesome
if numel(savevars) > 0
	[ST,I] = dbstack();
	save(ST(1).name,savevars{:});
end	
