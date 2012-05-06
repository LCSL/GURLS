function [X] = norm_l2(X,y,opt)
% norm_l2(X,Y,OPT)
%	Spheriphies the data according to the l2 norm.
% 
% INPUTS:
% -X: input data matrix
% -y: not used
% -OPT: not used
% 
% OUTPUT:
% -X: normalized input data matrix

for j = 1:size(X,1)
	X(j,:) = X(j,:)/(norm(X(j,:)) + eps);
end
fprintf('\tL2 Normalized\n');

