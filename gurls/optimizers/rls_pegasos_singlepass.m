function [cfr] = rls_pegasos_driver(X, bY, opt)
% rls_pegasos_singlepass(X,BY,OPT)
% utility function called by rls_pegasos
% computes a single pass for pegasos algorithm, performing the stochastic
% gradient descent over all training samples once.
%
% INPUTS:
% -X: input data matrix
% -BY: binary coded labels matrix
% -OPT: structure of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%       - epochs
%   fields that need to be added by hand
%       -Xte
%       -yte
%   fields with default values set through the defopt function:
%		- singlelambda
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: structure with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -W_sum: sum of the classifiers across iterations
% -t0: stepsize parameter 
% -count: number of iterations
% -acc_last: accuracy of the solution computed in the last iteration
% -acc_avg: average accuracy across iterations

lambda = opt.singlelambda(opt.paramsel.lambdas);

%% Inputs
[n,d] = size(X); 
[T] = size(bY,2);

%% Initialization
cfr = opt.cfr;

W = cfr.W;
W_sum = cfr.W_sum;
count = cfr.count;
t0 = cfr.t0;

%% Initialization
iter = 0;

seq = randperm(n); 
while iter < n,
    iter = iter + 1;
    idx = seq(iter);
    
    %% Stepsize
    %% Update Equations
    xt = X(idx,:);
    y_hat = (xt*W);
    r = bY(idx,:) - y_hat;
    %eta = 1.0/(lambda*(iter + t0));
    eta = 1.0/(lambda*(count + t0));
    W = (1 - lambda*eta)*W + eta*xt'*r;

    %% Projection onto the ball with radius sqrt(T/lambda)
    nW = norm(W,'fro');
    if nW > sqrt(T/lambda)
        W = (W/nW)*sqrt(T/lambda);
    end
    %% Averaging
    W_sum = W_sum + W;
    count = count + 1;

end
cfr.W = W;
cfr.W_last = W;
cfr.W_sum = W_sum;
cfr.count = count;
cfr.iter = iter;
cfr.C = [];
cfr.X = [];
end




