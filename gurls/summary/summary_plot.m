function [] = summary_plot(filestr, fields, nRuns, plotopt)
% Currently can't call it from  gurls.
% This opt is different from main opt.

for i = 1:numel(filestr)
	res{i} = summarize(filestr{i}, fields, nRuns{i});
end	
	% For each field, provide xlabel, ylabel and title of the plot
	
for f = 1:numel(fields)
	%figure(f);
	meanSideBySide = [];
	stdSideBySide = [];
	for i = 1:numel(filestr)
		meanSideBySide = [meanSideBySide res{i}{f}.perclassMean'];
		stdSideBySide  = [stdSideBySide res{i}{f}.perclassStd'];
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
    if isfield(plotopt,'legend')
        op.legend = plotopt.legend;
    else
        op.legend = filestr;
    end
	pretty_plot(meanSideBySide,stdSideBySide,op);
end

