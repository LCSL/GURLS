function [opt] = drls (X, y, opt)
% Currently does only linear kernel

lambda = opt.lambda;


fprintf('Solving dual RLS...\n');
n = size(opt.K,1);

K = opt.K + (n*lambda)*eye(n);

R = chol(K); 

opt.C = R\(R'\y);
