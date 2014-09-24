function [scores] = pred_primal(X,y, opt)
% pred_primal(opt)
% computes the predictions of the linear classifier stored in opt.rls.W 
% on the samples passed in the X matrix.
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%   -X: input data matrix
%   fields that need to be set through previous gurls tasks:
%		- rls.W (set by the rls_* routines)
% 
% OUTPUT:
% -scores: matrix of predicted labels
	scores = X*opt.rls.W;	
