function [] = pretty_plot(meanSideBySide,stdSideBySide,plotopt)
figure;
numbars = size(meanSideBySide,2);
handles.bars = bar(meanSideBySide);
hold on;
% This plots an error bar for each column
for i = 1:numbars
		x =get(get(handles.bars(i),'children'), 'xdata');
		x = mean(x([1 3],:));
		handles.errors(i) = errorbar(x, meanSideBySide(:,i), stdSideBySide(:,i), 'k', 'linestyle', 'none', 'linewidth', 1);
end
colormap 'summer';
if isfield(plotopt,'labels')
	set(gca,'XTickLabel', plotopt.labels);
	set(gca,'XLim',[0 numel(plotopt.labels)+2])
	set(gca,'XTick',[1:numel(plotopt.labels)])
	rotateticklabel(gca,45);
end	
if isfield(plotopt,'title')
	title(plotopt.title,'FontSize',16);
end	
if isfield(plotopt,'ylabel')
	ylabel(plotopt.ylabel,'FontSize',16);
end
if isfield(plotopt,'xlabel')
	xlabel(plotopt.xlabel,'FontSize',16);
end	
if isfield(plotopt,'legend')
	legend(plotopt.legend);
end	
