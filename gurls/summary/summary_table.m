function [meanSideBySide, stdSideBySide] = summary_table(filestr, fields, nRuns, plotopt)

for i = 1:numel(filestr)
	res{i} = summarize(filestr{i}, fields, nRuns{i});
end	
for f = 1:numel(fields)
	meanSideBySide = [];
	stdSideBySide = [];
	for i = 1:numel(filestr)
        meanSideBySide = [meanSideBySide res{i}{f}.perclassMean'];
		stdSideBySide =  [stdSideBySide res{i}{f}.perclassStd'];
	end	
	if isfield(plotopt,'titles')
		op.title = plotopt.titles{f};
	end
	if isfield(plotopt,'ylabels')
		op.ylabel = plotopt.ylabels{f};
	end
	if isfield(plotopt,'xlabels')
		op.xlabel = plotopt.xlabels{f};
	end
	if isfield(plotopt,'class_names')
		op.labels = plotopt.class_names;
	end
	op.legend = filestr;
    
    print_table(meanSideBySide,stdSideBySide,op);	
end	
	
