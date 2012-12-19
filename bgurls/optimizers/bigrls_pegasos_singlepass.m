function [cfr] = bigrls_pegasos_driver(X, bY, opt)
lambda = opt.singlelambda(opt.paramsel.lambdas);

%% Inputs
[n,d] = size(X); 
[T] = size(bY,2);

%% Initialization
cfr = opt.cfr;

W = cfr.W;
W_sum = cfr.W_sum;
count = cfr.count;
t0 = cfr.t0;

%% Initialization
iter = 0;

seq = randperm(n); 
while iter < n,
    iter = iter + 1;
    idx = seq(iter);
    
    %% Stepsize
    %% Update Equations
    xt = X(idx,:);
    y_hat = (xt*W);
    r = bY(idx,:) - y_hat;
    %eta = 1.0/(lambda*(iter + t0));
    eta = 1.0/(lambda*(count + t0));
    W = (1 - lambda*eta)*W + eta*xt'*r;

    %% Projection onto the ball with radius sqrt(T/lambda)
    nW = norm(W,'fro');
    if nW > sqrt(T/lambda)
        W = (W/nW)*sqrt(T/lambda);
    end
    %% Averaging
    W_sum = W_sum + W;
    count = count + 1;

%    %% Testing
%    if(mod(count,n) == 1)
%	
%    	fprintf('\n\tObjective : %f',obj_primal(W, X, bY, lambda));
%    	cfr.acc_last(end+1) = test_classifier (W,opt);
%    	fprintf('\n\tLast Acc: %f', cfr.acc_last(end));
%    	cfr.acc_avg(end+1) = test_classifier (W_sum/count,opt);
%    	fprintf('\n\tAvg Acc: %f\n', cfr.acc_avg(end));
%
%    end
 
    fprintf('%d\r',iter);
end
cfr.W = W;
cfr.W_last = W;
cfr.W_sum = W_sum;
cfr.count = count;
cfr.iter = iter;
cfr.C = [];
cfr.X = [];
end




