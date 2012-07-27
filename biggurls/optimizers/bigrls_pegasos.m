function [cfr] = bigrls_pegasos(X, bY, opt)

%	bigrls_pegasos(X,y,opt)
%	computes a classifier for the primal formulation of RLS.
%	The optimization is carried out using a stochastic gradient descent algorithm.
%	The regularization parameter is set to the one found in opt.paramsel (set by the bigparamsel_* routines).
%	In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
%	INPUT:
%		- X : input data bigarray
%		- Y : labels bigarray
%		- opt : struct witht he following fields:
%			- Fields set by other biggruls tasks:
%				* paramsel.lambdas (set by the bigparamsel_*) routines.
%			- Fields set through the bigdefopt function:
%				* singlelambda
%
%	OUTPUT: structure with the following fields:
%		- W : matrix of coefficient vectors of rls estimator for each class
%		- C : empty matrix
%		- X : empty matrix
				* epochs



	n = X.NumItems();
	d = X.Sizes();
	d = d{1};

	T = bY.Sizes();
	T = T{1};

opt.cfr.W = zeros(d,T);
opt.cfr.W_sum = zeros(d,T);
opt.cfr.count = 0;
opt.cfr.acc_last = [];
opt.cfr.acc_avg = [];


% Run mulitple epochs
for i = 1:opt.epochs,

		nBlocks = X.NumBlocks();
		blockOrder = randperm(nBlocks);
	
		for j = blockOrder,
	
			fprintf('Using block %d/%d\n',j,nBlocks);
			tX  = X.ReadBlock(j,'v');
		
			tbY = bY.ReadBlock(j,'v');
			blockBegin = (j-1)*X.BlockSize();
			blockSize = size(tX,2);
			tbY = bY(blockBegin+1: blockBegin+blockSize, :);
			
			if opt.cfr.count == 0
				opt.cfr.t0 = ceil(norm(tX(:,1))/sqrt(opt.singlelambda(opt.paramsel.lambdas)));
				fprintf('\n\tt0 is set to : %f\n', opt.cfr.t0);
			end
			% We expect to have a transposed bigarray;
			opt.cfr = bigrls_pegasos_singlepass(tX', tbY, opt);
		end
end	
cfr = opt.cfr;
cfr.W = opt.cfr.W_sum/opt.cfr.count;
