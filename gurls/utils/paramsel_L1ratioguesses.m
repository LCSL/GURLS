function [guesses] = paramsel_L1ratioguesses(y, opt, eig_list, tot1)
% paramsel_lambdaguesses(EIGVALS, R, N, OPT)
% Builds array of possible values for the regularization parameter, 
% generating a geometric series from the values in EIGVALS 
% Internal function, not to be called from gurls
% 
% INPUTS:
% -y: training label
% -eig_list: 1 x M double of largest eigenvalues of each kernel matrix
% -tot: number of candidate parameters
% -OPT: structure of options with the following fields with default values
%       set through the defopt function:
%		- mkl.smallnumber: use as lower limit of nonzero L1/L2 term
%
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUTS: 
% -GUESSES: 1 x tot array of candidate values for the L1 ratio 

eigK_app = max(eig_list);
M = length(eig_list);

lmax = (1/M) * norm(y)/(sqrt(eigK_app)) * 2;
lmin = opt.mkl.smallnumber;

% calculate range
power = linspace(log(lmin), 0, tot1);
guesses = lmax * exp(power);

guesses = [0, guesses];
end
