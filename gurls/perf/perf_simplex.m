function [p] = perf_simplex(X, y, opt)

%	perf_macroavg(X,y,opt)
% 	Computes the average classification accuracy per class.
%
%	NEEDS:
%		- opt.pred

if isstruct(opt.pred)
	opt.pred = opt.pred.means;
end	

if isfield (opt,'perf')
	p = opt.perf; % lets not overwrite existing performance measures.
		      % unless they have the same name
end


if(size(y,2)>1)
[dump,y_true] = max(y*opt.codes,[],2);

else 
y_true= y;

end
T=max(y_true);
y_pred = opt.pred;

%% Assumes single label prediction.
y_pred = opt.pred*opt.codes;

[dummy, predlab] = max(y_pred,[],2);
truelab 		 = y_true;
[MacroAvg, PerClass] = macroavg(truelab, predlab);

for t = 1:T,
	p.acc(t) = PerClass(t);
	p.forho(t) = p.acc(t);
	p.forplot(t) = p.acc(t);
end


function [MacroAverage, PerClass] = macroavg(TrueY, PredY)
% Computes average of performance for each class.

% Micro
% micro = mean(Classes == YTest);
% fprintf('Micro Avg: %2.4f\n',micro);

% Macro
nClasses = max(TrueY);
for i = 1:nClasses,
    acc(i) = sum((TrueY == i) & (PredY == i))/(sum(TrueY == i) + eps);
end
PerClass = acc;
MacroAverage = mean(acc);

%fprintf('Macro Avg: %2.4f\n',macro);
