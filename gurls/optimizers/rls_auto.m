function [cfr] = rls_auto(X, y, opt)
% Does automatic selection of primal/dual procedure.

% Very basic. Should think if we can make this more intelligent.

[n,d] = size(X);

if (n > d) % Do primal
	rls_primal(X, y, opt);
else % Do dual
	rls_dual(X, y, opt);
end
	
end
