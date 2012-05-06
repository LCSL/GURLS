function [scores] = pred_dual(X, y, opt)
% pred_primal(X,y,opt)
% computes the predictions of the classifier stored in opt.rls on 
% the samples passed in the X matrix.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: structure of options with the following fields (and subfields):
%   fields with default values set through the defopt function:
%		- kernel.type
%   fields that need to be set through previous gurls tasks:
%		- rls (set by the rls_* routines)
%       - predkernel.K (set by the predkernel_* routines), only if 
%         kernel.type is not 'linear'
% 
% OUTPUT:
% -scores: matrix of predicted labels

if strcmp(opt.kernel.type , 'linear')
	scores = pred_primal(X, y, opt);
else
	% Write the dual prediction code here.
	scores = opt.predkernel.K*opt.rls.C;
end
