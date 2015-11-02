function [fk] = predkernel_traintest_nystrom(X, y, opt)
% predkernel_traintest_nystrom(OPT)
% Computes the kernel matrix between the Nystrom subsampled training points and the
% test points. It can be used by pred_dual_nystrom.
%
% INPUTS:
% -OPT: structure with the following fields:
%   fields that need to be set through previous gurls tasks:
%       - X: input data matrix
%		- rls.X (only if opt.kernel.type is 'rbf' or 'chisquared' )
%		- testkernel (only if opt.kernel.type is 'load')
%       - paramsel.sigma (set by the paramsel_siglam* routines, is
%         required, only if opt.kernel.type is 'rbf')
%   fields with default values set through the defopt function:
%       - kernel.type
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with at least the field Knm containing the ntest x m kernel matrix

switch opt.kernel.type
    
	case {'rbf_nystrom'}
		fk.type = 'rbf_nystrom';
		if ~isprop(opt,'predkernel')
			opt.newprop('predkernel', struct());
			opt.predkernel.type = 'rbf_nystrom';
		end
        if ~isfield(opt.predkernel,'distance')
            opt.predkernel.distance = square_distance(X',opt.rls.X(opt.nystrom.sampledIdx,:)');
        end
		fk.distance = opt.predkernel.distance;
			D = -(opt.predkernel.distance);
			fk.Knm = exp(D/(opt.paramsel.sigma^2));
		if isfield(opt.rls,'L')
		    fk.Ktest = ones(size(X,1),1);
        end
        
    % TO DO: Support for kernels other than rbf_nystrom

    otherwise
        error('Kernel not supported yet');
end	
 
		
