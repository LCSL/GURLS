function [vout] = split_lap(X,y,opt);

[n,T] = size(y);

order = randperm(n);
last = floor(n*opt.hoproportion);
vout.va = order(1:last);
vout.tr = order(last+1:end);

l = floor((n-last)*opt.labelledproportion);

vout.l = vout.tr(1:l);
vout.u = vout.tr(l+1:end);


