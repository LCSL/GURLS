function [ gurls_chain ] = preproc_crossvalidate( X, y, gurls_chain, param_data)
%PREPROC_CROSSVALIDATE Crossvalidation for preprocessing.  This is also
%  probably pretty useful for other things since it is architected in a
%  general way.

% Inputs: 
%  X: The data. Nxd
%  y: The labels/output. Nx1
%  gurls_chain: a GURLS `opt` structure. This needs to be fairly specific:
%   (1) it needs to have a `train' process and a `test' process (processes
%   1 and 2, respectively)
%   (2) it needs to have a perf: structure
% param_data: Information for the parameter selection
%

% Right now I have random subsample cross validation implemented
n_trials = param_data.n_trials;
pct_train = param_data.pct_train;

% Current implementation only supports one parameter.
param_seq = regexp(param_data.param_name,'\.','split');
param_vals = param_data.param_vals;

n_params = length(param_data.param_vals);
perf_data = zeros(n_trials,n_params);


[N,d] = size(X);
for k_trial=1:n_trials
    ix = randperm(N)';
    ix_split = fix(N*pct_train);
    ix_tr = ix(1:ix_split);
    ix_te = ix((ix_split+1):end);
    fprintf('  [cross validation trial %d/%d]\n', k_trial,n_trials);
    for k_param=1:n_params
        % Set the parameter
        opt = deep_copy(gurls_chain);
        opt = set_param(opt,param_seq,param_vals(k_param));
        gurls(X(ix_tr,:),y(ix_tr),opt,1);
        gurls(X(ix_te,:),y(ix_te),opt,2);
        perf_data(k_trial,k_param) = opt.perf.(param_data.perf_field);
    end
end
mean(perf_data)
% Find maximum averaged over trials
switch param_data.max_or_min
    case 'max'
    [~,arg_opt] = max(mean(perf_data,1));
    case 'min'
    [~,arg_opt] = min(mean(perf_data,1));
end
best_values = param_vals(arg_opt);
gurls_chain = set_param(gurls_chain,param_seq,best_values);
return;
    function C = deep_copy(A)
        % This isn't great, but it's the only way to reliably deep copy
        % handle objects in MATLAB (that I know of).
        tmp = [tempname '.mat'];
        save(tmp,'A','-v7.3');
        load(tmp);
        C = A;
        delete(tmp);
    end

    function opt = set_param(opt,seq,val)
        if isempty(seq)
            opt = val;
        else
            opt = setfield(opt,seq{1},set_param(opt.(seq{1}),{seq{2:end}},val));
        end
    end
end

