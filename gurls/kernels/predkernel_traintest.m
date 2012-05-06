function [fk] = predkernel_traintest(X,y,opt)
% predkernel_traintest(X,Y,OPT)
% Computes the kernel matrix between the training points and the
% test points. It can be used by pred_dual.
%
% INPUTS:
% -X: input data matrix
% -Y: not used 
% -OPT: structure with the following fields:
%   fields that need to be set through previous gurls tasks:
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
% OUTPUT: ntestxntr matrix


switch opt.kernel.type
	case {'rbf'}
		fk.type = 'rbf';
		if ~isfield(opt,'predkernel')
			opt.predkernel.type = 'rbf';
		end
		if ~isfield(opt.predkernel,'distance')
			opt.predkernel.distance = distance(X',opt.rls.X');
			fk.distance = opt.predkernel.distance;
		end
		%D = -(opt.finalkernel.distance.^2);
		D = -(opt.predkernel.distance);
		fk.K = exp(D/(opt.paramsel.sigma^2));
	case {'load'}
		fk.type = 'load';
		load(opt.testkernel);
		fk.K = K_tetr;
	case {'chisquared'}
		for i = 1:size(X,1)
			for j = 1:size(opt.rls.X,1)
				fk.K(i,j) = sum(...
								( (X(i,:) - opt.rls.X(j,:)).^2 ) ./ ...
								( 0.5*(X(i,:) + opt.rls.X(j,:)) + eps));
			end
		end
end	
 
		
