function [scores] = pred_dual_nystrom(X,y, opt)
% pred_dual_nystrom(X,y,opt)
% computes predictions  of the input samples passed in the X matrix, 
% using the predictor stored in opt.rls.
%
% INPUTS:
%   -X: input data matrix
%   -y: labels matrix
%   -OPT: structure of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- rls (set by the rls_* routines)
%       - predkernel.Knm (set by the predkernel_*_nystrom routines)
% 
% OUTPUT:
% -scores: matrix of predicted labels

if isprop(opt, 'predkernel')
    scores = opt.predkernel.Knm * opt.rls.C;
else	
	error('Please provide a predkernel');
end

