function [] = plot_mkl_L1path(X, y, opt)
% plot_mkl_L1path(X, y, opt)
% take the output from paramsel_homkl, and plot 
% the L2 norm of dual parameters for each MKL kernel
% alone a grid of L1 penalty parameters, 
% with L2 parameter fixed at the optimal value (usually very close to 0)
% 
% INPUTS:
% -X, y: not used
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel (set by the paramsel_homkl routine)
%           - with following fields: path_mkl, guesses_mkl
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUTS: 
% -GUESSES: 1 x tot array of candidate values for the L1 ratio 
[guesses_L1, guesses_L2] = opt.paramsel.guesses_mkl{:};
M = size(opt.kernel.K)

tot1 = length(guesses_L1);

% solution path plot
clf
title('L1 regularization path for ||\alpha_m||_2 (m \in 1:M)')
xlabel('L1 penalty term^{1/3}') 
ylabel('||\alpha_m||_2') 

hold on
for par_id = 1:tot1
plot([guesses_L1(par_id), guesses_L1(par_id)].^(1/3),...
    [0, max(opt.paramsel.path_mkl(:))], ...
    'Color',[0.74,0.74,0.74]);
end
for m = 1:M
plot((guesses_L1).^(1/3), opt.paramsel.path_mkl(m, :));
end
hold off

end
