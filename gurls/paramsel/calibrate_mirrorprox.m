function [params] = paramsel_calibratesgd(X, y, opt)

%	paramsel_calibratesgd(X,y,opt)
%	Computes soem parameters needed by rls_smp (calibrate struct).
%
%	NEEDS:
%		- opt.smp_nestimates
%		- opt.smp_subsize

%% Subsample a set of 6000 examples;  We can change this later
%% Do 10 subsamples and ten estimates.
[n,d] = size(X);

T = size(y,2);

for i = 1:opt.smp_nestimates,

	idx = randsample(n, opt.smp_subsize);
	M = X(idx,:);

	if size(M,1) < size(M,2)
		K = M*M';
	else	
		K = (1/opt.smp_subsize)*(M'*M);
	end	
	%eigenmax(i) = max(eigs(K));
	rk(i) = rank(K);
	eigenmax(i) = normest(K);
%	eigenmax(i) = trace(K);
end	
params.c = 1/(1-sqrt(mean(rk(i))/opt.smp_subsize));
params.c = max(params.c,params.c*params.c);
params.eigenmax = params.c*mean(eigenmax);
%params.c = 1/(1-(3*sqrt(2)/sqrt(opt.smp_subsize)));
%params.eigenmax = params.c*mean(eigenmax);
params.gamma = 1/(sqrt(3)*params.eigenmax);

%s = zeros(opt.smp_nestimates,T);
%for k = 1:opt.smp_nestimates
%	
%	idx = randsample(n, opt.smp_subsize);
%	% This should work with both normal matrices
%	M = X(idx,:);
%	yk = y(idx,:);
%
%	w = M'*( (M*M')\yk );
%	% Maybe:
%	% w = M'*( ( M*M'+ eps*eye(size(M,1)) )\yk );
%
%	gk = -(1/opt.smp_subsize) * M' * (yk - M*w);
%	for p = 1:opt.smp_subsize
%		xp = M(p,:);
%		yp = yk(p,:)
%		vp = - xp'*(yp - xp*w);
%		sp = vp-gk;
%		np = sum(sp'.*sp,2); % same as: diag(sp'*sp);  
%		s(k,:) = s(k,:) + np;
%	end
%	s(k,:) = s(k,:)/opt.smp_subsize;
%end
%sbar = mean(s,2);
%shat = std(s,0,2);
