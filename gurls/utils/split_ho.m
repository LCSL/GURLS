function [vout] = split_ho(X,y,opt)
if size(X,1) > 0
	[n,d] = size(X);
else
	n = size(opt.kernel.K,1);
end	
order = randperm(n);
last = floor(n*opt.hoproportion);
vout.va = order(1:last);
vout.tr = order(last+1:end);

				
