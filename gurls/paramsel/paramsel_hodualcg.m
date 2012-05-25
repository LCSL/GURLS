function [vout] = paramsel_hodualcg(X,y,opt)

%	paramsel_hodualcg(X,y,opt)
%	Performs parameter selection when the dual formulation of conjugate gradient
% descent is used
%	The performance measure specified by opt.hoperf is maximized.
%
%	NEEDS:	
%		- opt.split
%		- opt.nholdouts
%		- opt.kernel.type
%		- opt.kernel.K
%		- opt.hoperf

if (isfield(opt, 'paramsel'))
	vout = opt.paramsel;
end

[n, T] = size(y); % number of samples, number of classes for multiclass

% Force opt.split to be a cell array, even if just using single holdout
if ~iscell(opt.split)
	s = opt.split;
	opt = rmfield(opt, 'split');
	opt.split{1} = s;
end
vout.ho_iter = zeros(1,opt.nholdouts);
K = opt.kernel.K;
maxiter = opt.gd.maxiter;
for j=1:opt.nholdouts
	tr = opt.split{j}.tr;
	va = opt.split{j}.va;

	  
	opt.kernel.K = (1/numel(tr))*K(tr,tr);
	opt.predkernel.K = K(va,tr);

	% Initialize
	opt.cg.a = zeros(numel(tr),T);
	opt.cg.r = y(tr,:);
	opt.cg.d = y(tr,:);
	opt.cg.t = opt.kernel.K*y(tr,:);

	
	ap = zeros(maxiter,T);
	
	% Iterate
	for i = 1:maxiter
		opt.cg = rls_dualcg_driver([],y(tr,:),opt);
		opt.rls.C = n*opt.cg.a;
		opt.pred = pred_dual([],[],opt);
		opt.perf = opt.hoperf([],y(va,:),opt);
	  	ap(i,:) = opt.perf.acc;
	  
	  	 % To do: Add early termination condition
	end

	[~, best_idx] = max(mean(ap,2));
	vout.ho_iter(j) = best_idx;
	vout.forho{j} =ap;
end

% Average over holdouts
vout.iter = mean(vout.ho_iter);
