function [fk] = predkernel_traintest(X, y, opt)
% predkernel_traintest(OPT)
% Computes the kernel matrix between the training points and the
% test points. It can be used by pred_dual.
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
% OUTPUT: struct with at least the field K containing the ntestxntr kernel matrix

switch opt.kernel.type
	case {'rbf'}
		fk.type = 'rbf';
		if ~isprop(opt,'predkernel')
			opt.newprop('predkernel', struct());
			opt.predkernel.type = 'rbf';
		end
        if ~isfield(opt.predkernel,'distance')
            opt.predkernel.distance = square_distance(X',opt.rls.X');
        end
		fk.distance = opt.predkernel.distance;
			%D = -(opt.finalkernel.distance.^2);
			D = -(opt.predkernel.distance);
			fk.K = exp(D/(opt.paramsel.sigma^2));
		if isfield(opt.rls,'L')
		    fk.Ktest = ones(size(X,1),1);
		end
	case {'load'}
		fk.type = 'load';
		load(opt.testkernel);
		fk.K = K_tetr;
    case {'linear'}
		fk.type = 'linear';
		fk.K = X*opt.rls.X';
	case {'chisquared'}
        for i = 1:size(X,1)
            for j = 1:size(opt.rls.X,1)
                fk.K(i,j) = sum(...
                    ( (X(i,:) - opt.rls.X(j,:)).^2 ) ./ ...
                    ( 0.5*(X(i,:) + opt.rls.X(j,:)) + eps));
            end
        end

	case {'datatype'}
		fk.type = 'datatype';
		fk.K = X;

    case {'quasiperiodic'}
		fk.type = 'quasiperiodic';
		if ~isprop(opt,'predkernel')
			opt.newprop('predkernel', struct());
			opt.predkernel.type = 'quasiperiodic';
		end
		if ~isfield(opt.predkernel,'distance')
			opt.predkernel.distance = square_distance(X',opt.rls.X');
		end
		fk.distance = opt.predkernel.distance;
		
		D = opt.predkernel.distance;
		fk.K = opt.paramsel.alpha*exp(-sin(((D.^(1/2)).*(pi/opt.period))).^2);
		fk.K = fk.K + (1-opt.paramsel.alpha)*exp(-D./(opt.paramsel.sigma^2));
		if isfield(opt.rls,'L')
		    fk.Ktest = ones(size(X,1),1);
        end
    
    otherwise
        kern = str2func(['kernel_', opt.kernelfun]);
        opt.kernel.sigma = opt.paramsel.sigma;
        opt.kernel.Y = opt.rls.X
        
        fk = kern(X, [], opt);
end	
 
		
