function [X] = norm_l2(X,y,opt)

%	norm_l2(X,y,opt)
%	Spheriphies the data according to the l2 norm.

for j = 1:size(X,1)
	X(j,:) = X(j,:)/(norm(X(j,:)) + eps);
end
fprintf('\tL2 Normalized\n');

