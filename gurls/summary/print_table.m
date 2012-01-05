function [] = print_table(meanSideBySide,stdSideBySide,op)

meanSideBySide =meanSideBySide;
stdSideBySide = stdSideBySide;
if isfield(op,'title')
	fprintf('\t%s\n',op.title);
end	
fprintf('\n');
fprintf('%8s\t','CLASS')
if isfield(op,'legend')
	for i = 1:numel(op.legend)
		fprintf('%s\t\t\t',op.legend{i});
	end
end
fprintf('\n\n\t\t');
for i = 1:size(meanSideBySide,2)
	fprintf('%8s\t%8s\t','MEAN','STD');
end	
fprintf('\n\n');
for i = 1:size(meanSideBySide,1) 		% Loop over classes
	if isfield(op,'labels')
		fprintf('%8s\t',op.labels{i});
	end	
	for j = 1:size(meanSideBySide,2)	% Loop over files
		fprintf('%f\t%f\t',meanSideBySide(i,j),stdSideBySide(i,j));
	end
	fprintf('\n');
end	
fprintf('\n\n');
