function [scores] = pred_dual(X,y, opt)
% pred_primal(X,y,opt)
% computes the predictions of the classifier stored in opt.rls on 
% the samples passed in the X matrix.
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%   -X: input data matrix
%   -y: labels matrix
%   fields with default values set through the defopt function:
%		- kernel.type
%   fields that need to be set through previous gurls tasks:
%		- rls (set by the rls_* routines)
%       - predkernel.K (set by the predkernel_* routines), only if 
%         kernel.type is not 'linear'
% 
% OUTPUT:
% -scores: matrix of predicted labels

if ~strcmp(opt.kernel.type, 'linear') && isprop(opt, 'predkernel')
	scores = opt.predkernel.K*opt.rls.C;
elseif strcmp(opt.kernel.type, 'linear')
	scores = pred_primal(X, y, opt);
else	
	error('Please provide a predkernel');
end
