function [params] = bigparamsel_calibratesgd(X, y, opt)

%   bigparamsel_calibratesgd(X, y, opt)
%	Performs parameter selection when one wants to solve the problem using rls_pegasos.
%
%	INPUT:
%		- X : input data bigarray
%		- Y : labels bigarray
%		- OPT : struct witht he following fields:
%			- Fields set through the bigdefopt function:
%				* nlambda
%				* hoperf
%               * subsize
%               * calibfile
%               * hoperf
%               * singlelambda
%			- Fields that need to be set by hand:
%				* opt.files.Xva_filename 	: Validation data bigarray
%				* opt.files.yva_filename 	: Validation labels bigarray
%
%   For more information on standard OPT fields
%   see also defopt
% 
%   OUTPUTS: structure with the following fields:
%       - lambdas:  selected value for the regularization parameter
%       - W:        rls coefficient vector

n_estimates = 1;
	n = X.NumItems();

for i = 1:n_estimates,

	sub_size = opt.subsize;
	idx = randsample(n, sub_size);

	M = X(idx,:);


	if ~exist([opt.calibfile '.mat'],'file')
		fprintf('\n\tCalibrating...');
		%% Step 1 : Hold out parameter selection in the dual
		name = [opt.calibfile];
		tmp.hoperf = opt.hoperf;
		tmp = defopt(name);
		tmp.seq = {'split:ho','kernel:linear','paramsel:hodual','rls:dual'};
		tmp.process{1} = [2,2,2,2];
		tmp.singlelambda = opt.singlelambda;

		gurls(M,y(idx,:),tmp,1);
	end
	fprintf('\n\tLoading existing calibration');
	load([opt.calibfile '.mat']);
	lambdas(i) = opt.singlelambda(opt.paramsel.lambdas);

	% Add rescaling
end
params.lambdas = mean(lambdas);
params.W = opt.rls.W;

%% Step 2 : Calibrate stepsize parameter t0
% This is done inside the pegasos code.

%% Step 3 : Calibrate the Stepsize parameter 1/lambda*(t + t0)
%% Which value of t0 gets maximum reduction in the objective function?


%% Step 3 : Calibrate the Stepsize parameter 1/lambda*(t + t0)
%% Which value of t0 gets maximum reduction in the objective function?


