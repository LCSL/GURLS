function [kernel] = kernel_quasiperiodic(X,y,opt)
% 	kernel_quasiperiodic(X,y,opt)
%	Computes the kernel matrix for a univariate quasiperiodic kernel, 
%   convex sum of periodic and gaussian kernel.
%	INPUTS:
%		-X: input data matrix
%		-y: not used 
%		-OPT: struct with the following options:
%       - paramsel : struct containing the following fields (computed by paramsel_*).
%		- sigma : width of the gaussian kernel.
%	
%	OUTPUT: struct with the following fields:
%		-type: 'quasiperiodic'
%		-K: kernel matrix


if ~isfield(opt,'kernel')
	opt.kernel.type = 'quasiperiodic';
end


if ~isfield(opt.kernel,'distance')
	opt.kernel.distance = square_distance(X',X');
    n = length(X);
end	

kernel = opt.kernel;

D = opt.kernel.distance;
K = opt.paramsel.alpha*exp(-sin(((opt.kernel.distance.^(1/2)).*(pi/opt.period))).^2);
K = K+ (1-opt.paramsel.alpha)*exp(-D./(opt.paramsel.sigma^2));

kernel.type = 'quasiperiodic';
kernel.K = K;
