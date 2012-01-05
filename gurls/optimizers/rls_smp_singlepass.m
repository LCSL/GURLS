function [cfr] = rls_smp_singlepass(X,y,opt)

cfr = opt.cfr;
w = cfr.W;
r = cfr.R;

wsum = cfr.W_sum;
count = cfr.count;

gamma = opt.calibrate.gamma

lambda = opt.singlelambda(opt.paramsel.lambdas);

[n,d] = size(X); 
T = size(y,2);
seq = randperm(n);
for i = 1:n

	count = count + 1;

	x = X(seq(i),:);
	yy = y(seq(i),:);

	w = r + gamma * (x') * (yy - x * r);
	%proj
	nW = norm(w,'fro');
    if nW > sqrt(T/lambda)
        w = (w/nW)*sqrt(T/lambda);
    end

	r = r + gamma * (x') * (yy - x * w);
	%proj
	
	%wsum = wsum + gamma*w
	%wavg = wsum/(gamma*count);

	wsum = wsum + w;

	if(mod(count,n) == 1)
    	fprintf('\n\tObjective : %f',obj_primal(w, X, y, 0));
		if isfield(opt,'Xte') && isfield(opt,'yte')
	    	cfr.acc_last(end+1) = test_classifier (w,opt);
   		 	fprintf('\n\tLast Acc: %f', cfr.acc_last(end));
    		cfr.acc_avg(end+1) = test_classifier (wsum/count,opt);
    		fprintf('\n\tAvg Acc: %f\n', cfr.acc_avg(end));
		end	
    end
 
    fprintf('%d\r',i);
end

cfr.W = w;
cfr.W_sum = wsum;
cfr.R = r;
cfr.count = count;
