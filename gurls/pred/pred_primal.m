function [scores] = pred_primal(X, y, opt)
% pred_primal(X,y,opt)
% computes the predictions of the linear classifier stored in opt.rls.W 
% on the samples passed in the X matrix.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: structure of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- rls.W (set by the rls_* routines)
% 
% OUTPUT:
% -scores: matrix of predicted labels

	scores = X*opt.rls.W;	
