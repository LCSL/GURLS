function [cfr] =rls_svmcons(X,Y,opt)
% y sent is encoded output 
acc_avg 		= 0;
acc_last 		= 0;
epochs 			= opt.epochs;
maxIter 		= opt.maxIter;
sampling 		= opt.sampling;
batchsize 		= opt.batchsize;
lambda 			= opt.singlelambda(opt.paramsel.lambdas);
t				= opt.t;

%% Inputs

[n,d]  			= size(X);
coding 			= opt.codes';
W      			= zeros(t-1,d);


%% Setup
maxIter 		= max(maxIter,opt.epochs*n);
iter 			= 0;
W_sum 			= W;
count 			= 0;
t0				= 0;

while iter < maxIter

 
    iter = iter + 1;
    %% Epochs
    switch sampling
        case {'epoch'}
            %  Shuffle once at the beginning of the epoch
            if (mod(iter,n) == 1), order = randperm(n); 
      
	    %fprintf('\nEpoch: %d', floor(iter/n)); 
		end
            idx = order(mod(iter,n) + 1);
        case {'random'}
            %% Random sampling (No Epochs)
            idx = randsample(n,batchsize);
        otherwise
            error(['Unknown sampling type: ', options.sampling]);
    end
    
    xt = X(idx,:);
    y_hat = (W*xt'); 
    
    eta = 1.0/(lambda*((iter)));

    projt = coding*y_hat+1/(t-1);

	coding(Y(idx),:) = 0;
	r = sum(coding(projt > 0,:),1);
	coding = opt.codes';

	if coding(Y(idx),:)*y_hat < 1
	    W = (1 - lambda*eta)*W - eta*((r-coding(Y(idx),:))'*xt);
	else
	    W = (1 - lambda*eta)*W - eta*(r'*xt);
	end
    nW = norm(W,'fro');
    if nW > sqrt(1/lambda)
        W = W/(nW*sqrt(lambda));
    end
    if iter > 0, % Start averaging
       W_sum = W_sum + W;
       count = count + 1;
    end
%    if(mod(iter,n) == 1)
%    	fprintf('\n\tNorm of W : %f, %f',nW, sqrt(1/lambda));
%    	acc_last(end+1) = test_classifier (W,opt);
%    	fprintf('\n\tLast Acc: %f', acc_last(end));
%    	acc_avg(end+1) = test_classifier (W_sum/count,opt);
%    	fprintf('\n\tAvg Acc: %f\n', acc_avg(end));
%        c=X*(W_sum/count)';
%       
%    end
end

cfr.W = (W_sum'/count);
cfr.W_avg = (W_sum/count)';
cfr.C = [];
cfr.X = [];
cfr.acc_last = acc_last;
cfr.acc_avg = acc_avg;


end

function [acc] = test_classifier (W, opt)
	opt.rls.W = W';
	opt.pred = pred_primal(opt.Xte, opt.yte, opt);
	opt.perf   = perf_simplex(opt.Xte, opt.yte, opt);
	acc = mean([opt.perf.acc]);
end
