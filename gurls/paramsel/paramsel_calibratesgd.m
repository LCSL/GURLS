function [params] = paramsel_calibratesgd(X,y, opt)
% paramsel_calibratesgd(X,y, OPT)
% Performs parameter selection when one wants to solve the problem using rls_pegasos.
%
% INPUTS:
% -OPT: structure of options with the following fields:
%   with default values  set through the defopt function:
%		- subsize
%		- calibfile
%		- hoperf
%		- singlelambda
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUTS: structure with the following fields:
% - lambda: selected value for the regularization parameter
% - W: rls coefficient vector
%

n_estimates = 1;
[n,d] = size(X);



for i = 1:n_estimates,

	sub_size = opt.subsize;
    if n<sub_size;
        error('GURLS usage error: the option subsize of the option list must be smaller than the number of training samples!!');
    end
	idx = randsample(n, sub_size);

	M = X(idx,:);


	if ~exist([opt.calibfile '.mat'],'file')
        if opt.verbose
    		fprintf('\n\tCalibrating...');
        end
		%% Step 1 : Hold out parameter selection in the dual
		name = opt.calibfile;
		tmp = defopt(name);
		tmp.hoperf = opt.hoperf;
		tmp.seq = {'split:ho','kernel:linear','paramsel:hodual','rls:dual'};
		tmp.process{1} = [2,2,2,2];
		tmp.singlelambda = opt.singlelambda;

		gurls(M,y(idx,:),tmp,1);
    end
    if opt.verbose
    	fprintf('\n\tLoading existing calibration');
    end
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


