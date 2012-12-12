function [scores] = pred_randfeats(X, y, opt)
% pred_primal(X,y,opt)
% computes the predictions of the linear classifier stored in opt.rls.W, 
% and obtained the Random Features approach proposed in: 
%   Ali Rahimi, Ben Recht;
%   Random Features for Large-Scale Kernel Machines;
%   in Neural Information Processing Systems (NIPS) 2007.
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

    G = rp_apply_real(X',opt.kernel.W)';
	scores = G*opt.rls.W;	
