function [scores] = pred_primal(X,y,opt)

%	bigpred_primal(X,y,opt)
%	computes the predictions of the linear classifier stored in opt.rls.W (generated
%	by some bigrls_* method) on the points passed in the X matrix.
%
%	NEEDS:
%		- opt.rls.W
%		- opt.files.pred_filename
%		- opt.blocksize

		scores = bigarray_mat(opt.files.pred_filename);
		scores.Clear();
		scores.Init(opt.blocksize);
		for i = 1:X.NumBlocks();
			M = X.ReadBlock(i);
			scores.Append(opt.rls.W'*M);
		end	
		scores.Flush();
end		
