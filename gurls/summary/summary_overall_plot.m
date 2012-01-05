function [] = summary_overall_plot(filestr, fields, nRuns, plotopt)
% Currently can't call it from  gurls.
% This opt is different from main opt.
for i = 1:numel(filestr)
	res{i} = summarize(filestr{i}, fields, nRuns{i});
end	
for f = 1:numel(fields)
	figure;
	meanSideBySide = [];
	stdSideBySide = [];
	for i = 1:numel(filestr)
		meanSideBySide = [meanSideBySide res{i}{f}.overallMean];
		stdSideBySide =  [stdSideBySide  res{i}{f}.overallStd];
	end	
	bar(meanSideBySide);
	hold on;
	h = errorbar(meanSideBySide,stdSideBySide,'xr');
	colormap 'summer';
    
    if isfield(plotopt,'legend')
        legend_str = plotopt.legend;
    else
        legend_str = filestr;
    end
    
	set(gca,'XLim',[0 numel(legend_str)+2])
	set(gca,'XTick',[1:numel(legend_str)])
	set(gca,'XTickLabel', legend_str);
	rotateticklabel(gca,45);
	set(gca,'FontSize',16);

	if isfield(plotopt,'titles')
		title(plotopt.titles{f});
	end
	if isfield(plotopt,'ylabels')
		ylabel(plotopt.ylabels{f});
	end
	if isfield(plotopt,'xlabels')
		xlabel(plotopt.xlabels{f});
	end
end
