function [cfr] = rls_dual_mkl (X, y, opt)

% rls_dual_mkl(X, y, opt)
% computes a predictor based on MKL framework with elastic net regularization.
%
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -OPT: struct of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- kernel.K_mkl (set by the kernel_mkl routine)
%		- paramsel.par_mkl (set by the paramsel_homkl routine)
%   fields that need to be set by hand:
%		- rls.C_init (initial value for C_init)
%		- mkl.par_mkl (user specified L1/L2 parameter)
%   For more information on standard OPT fields
%   see also defopt
%
% OUTPUT: struct 'crf'(opt.rls) with the following fields:
% -W: empty matrix
% -C_mkl: n x n x M matrix of coefficient vectors of MKL rls estimator
% -X: training samples used by the routine

if isprop(opt,'rls')
	cfr = opt.rls; % lets not overwrite existing parameters.
			      		 % unless they have the same name
else
    opt.newprop('rls', struct());
end

% prepare parameters
verbose = opt.mkl.verbose.rls;

eig_app = max(opt.kernel.eig_mkl);
verbose = opt.mkl.verbose.rls;
iter_max = opt.mkl.iter_max.rls;
crit = opt.mkl.crit.rls;

if isfield(opt.mkl, 'par_mkl')
    [par_L1, par_L2] = opt.mkl.par_mkl{:};
else
    [par_L1, par_L2] = opt.paramsel.par_mkl{:};
end

if isfield(cfr, 'C_init')
    C_init = opt.rls.C_init;
else
    C_init = [];    
end

[cfr.C_mkl, cfr.epath_mkl] = rls_dual_mkl_pfbs(...
    opt.kernel.K_mkl, y, par_L1, par_L2, ...
    C_init, eig_app, iter_max, crit, verbose);

cfr.W = [];
cfr.X = X;
end
