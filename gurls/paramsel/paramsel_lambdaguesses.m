function [guesses] = paramsel_lambdaguesses(eigvals, r, n, opt)
% Internal function; 
% Not to be called from gurls;

% eigvals = a vector containing the eigenvalues of X'X or XX'
% r = rank; usually send min(n,d)

% ensuring eigenvalues are sorted
eigvals = sort(eigvals,'descend');

% maximum eigenvalue
lmax = eigvals(1);

% just in case, when r = min(n,d) and r x r has some zero eigenvalues
% we take a max; 200*sqrt(eps) is the legacy number used in the previous
% code, so i am just continuing it.
lmin = max(min(opt.smallnumber, eigvals(r)), 200*sqrt(eps));

powers = linspace(0,1,opt.nlambda);

guesses = lmin.*(lmax/lmin).^(powers);

guesses = guesses/n;

end
