function [guesses] = paramsel_lambdaguesses(eigvals, r, n, opt)
% paramsel_lambdaguesses(EIGVALS, R, N, OPT)
% Builds array of possible values for the regularization parameter, 
% generating a geometric series from the values in EIGVALS 
% Internal function, not to be called from gurls
% 
% INPUTS:
% -EIGVALS: a vector containing the eigenvalues of X'X or XX' where X is
% the input data matrix
% -R: rank; usually send min(n,d)
% -N: number of samples
% -OPT: structure of options with the following fields with default values
%       set through the defopt function:
%		- nlambda
%		- smallnumber
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUTS: 
% -GUESSES: opt.nlambdaX1 array of values for the regularization parameter

% ensuring eigenvalues are sorted
eigvals = sort(eigvals,'descend');

% maximum eigenvalue
lmax = eigvals(1);

% just in case, when r = min(n,d) and r x r has some zero eigenvalues
% we take a max; 200*sqrt(eps) is the legacy number used in the previous
% code, so i am just continuing it.

lmin = max(min(lmax*opt.smallnumber, eigvals(r)), 200*sqrt(eps));

powers = linspace(0,1,opt.nlambda);
guesses = lmin.*(lmax/lmin).^(powers);
guesses = guesses/n;

end
