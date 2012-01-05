function [] = summary_overall_table(filestr,fields,nRuns,op)
for i = 1:numel(filestr)
	res{i} = summarize(filestr{i}, fields, nRuns{i});
end	
for f = 1:numel(fields)
	meanSideBySide = [];
	stdSideBySide = [];
	for i = 1:numel(filestr)
		meanSideBySide = [meanSideBySide res{i}{f}.overallMean];
		stdSideBySide =  [stdSideBySide  res{i}{f}.overallStd];
	end	

	%meanSideBySide = 100*meanSideBySide;
	%stdSideBySide = 100*stdSideBySide;

	if isfield(op,'titles')
		fprintf('\t%s\n',op.titles{f});
	end	
	fprintf('\n');
	for i = 1:numel(filestr)
		fprintf('%8s\t',filestr{i});
	end
	fprintf('\n\n');
	for i = 1:size(meanSideBySide,2)
		fprintf('%s\t%s\t','MEAN','STD');
	end	
	fprintf('\n\n');
	for j = 1:size(meanSideBySide,2)	% Loop over files
		fprintf('%5.3f\t%5.3f\t',meanSideBySide(j),stdSideBySide(j));
	end
	fprintf('\n\n');
end	

