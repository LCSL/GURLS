function [] = summary_overall_plot(filestr, fields, nRuns, plotopt, fileroot)
% (?) Currently can't call it from gurls. This opt is different from main opt.
%
% INPUT:
% - plotopt: stores xlabel, ylabel and title of the plot for each field.
% - fileroot: location of files pointed by filestr. Default is '.'

if ~exist('fileroot', 'var')
    fileroot = '.';
end

for indFile = 1:numel(filestr)
	res{indFile} = summarize(fullfile(fileroot, filestr{indFile}), fields, nRuns(indFile));
end

for f = 1:numel(fields)
	
	meanSideBySide = [];
	stdSideBySide = [];
	for i = 1:numel(filestr)
		meanSideBySide = [meanSideBySide res{i}{f}.overallMean];
		stdSideBySide =  [stdSideBySide  res{i}{f}.overallStd];
	end	
    
    figure;
	h1 = bar(meanSideBySide'); hold on;
    set(h1, 'FaceColor', summer(1), 'EdgeColor', 'k');
	h = errorbar(meanSideBySide, stdSideBySide, 'xk');
	% colormap 'summer'; % axis tight; fn = 12;
    
    if isfield(plotopt,'legend')
        legend_str = plotopt.legend;
    else
        legend_str = filestr;
    end
    
    fn = 12; % fontsize
	set(gca,'XLim', [0 numel(legend_str)+2])
	set(gca,'XTick', (1:numel(legend_str)))
	set(gca,'XTickLabel', legend_str);
	rotateticklabel(gca, 45);
	set(gca, 'FontSize', fn);
    
	if isfield(plotopt,'titles')
		title(plotopt.titles{f}, 'FontSize', fn);
	end
	if isfield(plotopt,'ylabels')
		ylabel(plotopt.ylabels{f}, 'FontSize', fn);
	end
	if isfield(plotopt,'xlabels')
		xlabel(plotopt.xlabels{f}, 'FontSize', fn);
	end
end
