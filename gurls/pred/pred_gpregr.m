function [pred] = pred_gpregr(X,y, opt)
% pred_gpregr(opt)
% Computes the predictive distribution (mean and variance) for a gaussian
% process
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%   -y: labels matrix
% 
% OUTPUT:  struct with the following fields:
%  -means
%  -vars

pred.means = opt.predkernel.K*opt.rls.alpha;
n = size(y,1);
pred.vars = zeros(n,1);
for i = 1:n;
    v = opt.rls.L'\opt.predkernel.K(i,:)';
    pred.vars(i) = v'*v; 
end
pred.vars = opt.predkernel.Ktest - pred.vars;
