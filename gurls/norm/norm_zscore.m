function [X] = norm_zscore(X,y, opt)
% norm_zscore(OPT)
% Normalizes the data, centering them and rescaling them so that each dimension has std. deviation 1.
% Then saves stats in a file with name root specified in opt.name 
%
% INPUTS:
% -OPT: structure of options with the following fields (and subfields):
%		- opt.name
%       - opt.X
% 
% OUTPUT:
% -X: normalized input data matrix

savevars = {'meanX','stdX'};

[n,d] = size(X);

meanX = mean(X);
stdX = std(X) + eps;

X = X - repmat(meanX, n, 1);
X = X./repmat(stdX, n, 1);

if numel(savevars) > 0
	[ST,I] = dbstack();
	save([opt.name '_' ST(1).name],savevars{:}, '-v7');
end
