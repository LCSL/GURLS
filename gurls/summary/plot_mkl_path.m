function [] = plot_mkl_path(X, y, opt, outcome)
% plot_mkl_L1path(X, y, opt)
% take the output from paramsel_homkl, and plot
% the L2 norm of dual parameters for each MKL kernel
% alone a grid of L1 penalty parameters,
% with L2 parameter fixed at the optimal value (usually very close to 0)
%
% INPUTS:
% -X, y: not used
% -outcome: 'norm' for kernel estimate norm, 'perf' for performance
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel (set by the paramsel_homkl routine)
%           - paramsel.norm_path (nL1 x M x nfold)
%           - paramsel.perf_path (nL1 x M x nfold)
% 			- paramsel.guesses_mkl (cell array of candidate L1/L2 parameters)
%   For more information on standard OPT fields
%   see also defopt
%
% OUTPUTS:
% -GUESSES: 1 x tot array of candidate values for the L1 ratio

guesses_L1 = opt.paramsel.guesses_mkl{:};
M = size(opt.kernel.K_mkl, 3);
tot1 = length(guesses_L1);

if strcmp(outcome, 'norm')
    norm_path = opt.paramsel.norm_path(:,:,1)';
    i1_opt = find(opt.paramsel.guesses_mkl{1} == ...
        opt.paramsel.par_mkl{1});
    
    % solution path plot
    clf
    title('L1 regularization path for ||\alpha_m||_2 (m \in 1:M)')
    xlabel('L1 penalty term^{1/3}')
    ylabel('||\alpha_m||_2')
    
    hold on
    for par_id = 1:tot1
        % grid of L1 parameters used
        plot([guesses_L1(par_id), guesses_L1(par_id)].^(1/3),...
            [0, max(norm_path(:))], 'Color',[0.74,0.74,0.74]);
    end
    %{
    % highlight selected L1 parameter
    plot([guesses_L1(i1_opt), guesses_L1(i1_opt)].^(1/3),...
        [0, max(norm_path(:, i1_opt))], ...
        'Color',[1,0,0], 'LineWidth', 2);
    %}
    for m = 1:M
        % solution path
        plot((guesses_L1).^(1/3), norm_path(m, :));
    end
    hold off
elseif  strcmp(outcome, 'perf')
    perf_name = fieldnames(opt.perf);
    perf_path = mean(opt.paramsel.perf_path, 2);
    if strcmp(perf_name{1}, 'rmsestd')
        perf_path = -perf_path;
    end
    plot(guesses_L1.^(1/3), perf_path);
    title(sprintf('L1 regularization path for performance (%s)', ...
        perf_name{1}));
    xlabel('L1 penalty term^{1/3}')
    ylabel(perf_name{1})

end
end
