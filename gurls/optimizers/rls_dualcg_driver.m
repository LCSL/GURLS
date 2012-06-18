function cg = rls_dualcg_driver(X,y,opt)

% utility function called by rls_dualcg
% computes a single step for conjucate gradient descent algorithm
% NEEDS:
%	- opt.cg (as initialized by rls_dualcg)
%	- opt.kernel.K
%	- opt.kernel.type

T = size(y,2);
n = size(opt.kernel.K,1);
K = (1/n) * opt.kernel.K;
for class = 1:T
	a = opt.cg.a(:,class);
	r = opt.cg.r(:,class);
	d = opt.cg.d(:,class);
	t = opt.cg.t(:,class);

	nt = (1/n)*(t'*K*t);
	t = t/nt;
	d = d/nt;
	gamma = (1/n) * y(:,class)'*K*t;
	a = a + gamma * d;
	r = r - gamma * t;
	d = r -(1/n) * d * (t'*K*K*r);
	t = K * d;

	cg.a(:,class) = a;
	cg.r(:,class) = r;
	cg.d(:,class) = d;
	cg.t(:,class) = t;
end	


