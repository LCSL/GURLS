function [scores] = pred_primal(X,y,opt)

% bigpred_primal(X,y,opt)
% computes the predictions of the linear classifier stored in opt.rls.W 
% on the samples passed in the X matrix.
%
% INPUTS:
% -X: input data bigarray
% -y: labels bigarray
% -OPT: structure of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%		- rls.W (set by the rls_* routines)
%	- Fields that need to be set by hand:
%
%		* opt.files.pred_filename 	: prediction bigarray filename

% 
% OUTPUT:
% -scores: bigarray with predicted scores.


		scores = bigarray_mat(opt.files.pred_filename);
		scores.Clear();
		scores.Init(opt.blocksize);
		for i = 1:X.NumBlocks();
			M = X.ReadBlock(i);
			scores.Append(opt.rls.W'*M);
		end	
		scores.Flush();
end		
