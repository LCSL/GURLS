function [vout] = paramsel_hodualgd(X, y, opt)

%       paramsel_hodualgd(X,y,opt)
%       Performs parameter selection when the dual formulation of gradient
% descent is used
%       The performance measure specified by opt.hoperf is maximized.
%
%       NEEDS:
%   - opt.gd.method (0 for standard gd, 1 for accelerated)
%   - opt.gd.maxiter
%               - opt.split
%               - opt.nholdouts
%               - opt.kernel.type
%               - opt.kernel.K
%               - opt.hoperf

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

% Get step size
if (opt.gd.method == 0)
       max_eig = eigs(double(opt.kernel.K),1);
       eta = opt.gd.eta_numerator/max_eig;
elseif (opt.gd.method == 1)
   if (strcmp(opt.kernel.type,'rbf'))
	   eta=1/n;
   else
	   max_eig = eigs(double(opt.kernel.K),1);
	   eta=1/(max_eig*n);%good for linear
   end
       %eta = opt.gd.eta_numerator/max_eig;

   nu = opt.gd.nu;
else
       error('invalid opt.gd.method')
end
vout.eta = eta;

vout.ho_iter = zeros(1,opt.nholdouts);
K = opt.kernel.K;
for j=1:opt.nholdouts
       tr = opt.split{j}.tr;
       va = opt.split{j}.va;


       opt.kernel.K = K(tr,tr);
       opt.predkernel.K = K(va,tr);

       % Initialize
       coeffs = zeros(length(tr),T);
       if (opt.gd.method == 1)
               opt.gd.alpha1 = zeros(length(tr),T);
       end
       opt.gd.c = coeffs;
       opt.paramsel.eta = eta;

       ap = zeros(opt.gd.maxiter,T);

       % Iterate
       for i = 1:opt.gd.maxiter
               opt.gd.iter = i;
               opt.gd = rls_dualgd_driver([],y(tr,:),opt);
               opt.rls.C = opt.gd.c;
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
