function vout = paramsel_bfprimal(X,y,opt)
% paramsel_bfprimal(X,y, OPT)
% Performs parameter selection when the primal formulation of RLS is used.
% This method uses the hold-out cross validation approach in a brute force way
% i.e. the RLS problem is solved from scratch for each value of the regularizer.
%
% INPUTS:
% -OPT: structure of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- split (set by the split_* routine)
%   fields with default values set through the defopt function:
%		- hoperf
%       	- nholdouts
%  fields that need to be set by hand
%		- opt.paramsel.guesses : values for lambda one wishes to try
%		- opt.paramsel.optimizer : algorithm to solve the RLS problem
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -lambdas_round: cell array (opt.nholdoutsX1). For each split a cell contains the 
%       values of the regularization parameter lambda minimizing the 
%       validation error for each class.
% -forho: cell array (opt.nholdoutsX1). For each split a cell contains a matrix 
%       with the validation error for each lambda guess and for each class
% -guesses: cell array (opt.nholdoutsX1). For each split a cell contains an 
%       array of guesses for the regularization parameter lambda
% -lambdas: mean of the optimal lambdas across splits

if isprop(opt,'paramsel')
	vout = opt.paramsel; % lets not overwrite existing parameters.
			      		 % unless they have the same name
else
    opt.newprop('paramsel', struct());
end

[n,T] = size(y);

for nh = 1:opt.nholdouts
	if strcmp(class(opt.split),'cell')
		tr = opt.split{nh}.tr;
		va = opt.split{nh}.va;
	else	
		tr = opt.split.tr;
		va = opt.split.va;
	end	

	for i = 1:numel(opt.paramsel.guesses)
		opt.paramsel.lambdas = opt.paramsel.guesses(i);
		opt.newprop('rls', opt.paramsel.optimizer(X(tr,:),y(tr,:),opt));
		opt.newprop('pred', pred_primal(X(va,:),y(va,:),opt));
		opt.newprop('perf', opt.hoperf(X(va,:),y(va,:),opt));
		for t = 1:T
			ap(i,t) = opt.perf.forho(t);
		end	
	end	
	[dummy,idx] = max(ap,[],1);	
	vout.lambdas_round{nh} = opt.paramsel.guesses(idx);
	vout.forho{nh} = ap;

end
if numel(vout.lambdas_round) > 1
	lambdas = cell2mat(vout.lambdas_round');
	vout.lambdas = mean(lambdas);
else
	vout.lambdas = vout.lambdas_round{1};
end
