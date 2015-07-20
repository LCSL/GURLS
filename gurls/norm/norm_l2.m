function [X] = norm_l2(X,y, opt)
% norm_l2(OPT)
%	Spheriphies the data according to the l2 norm.
% 
% INPUTS:
% -OPT: must contain the field X
% 
% OUTPUT:
% -X: normalized input data matrix

for j = 1:size(X,1)
	X(j,:) = X(j,:)/(norm(X(j,:)) + eps);
end
if opt.verbose
    fprintf('\tL2 Normalized\n');
end

