function pretty_plot(meanSideBySide, stdSideBySide, plotopt)

fn = 12; % fontsize
lw = 1.5; % linewidth

figure;
[numgroups, numbars] = size(meanSideBySide);
handles.bars = bar(meanSideBySide);
hold on;

% Plot an error bar for each column
for b = 1:numbars
    try % Compatibility: older version that does not seem to work
        x = get(get(handles.bars(b), 'children'), 'xdata');
        x = mean(x([1 3],:));
    catch
        % x = handles.bars(b).('XData');
        groupwidth = min(handles.bars(b).('BarWidth'), numbars/(numbars + 1.5));
        x = (1:numgroups) - groupwidth/2 + (2*b-1) * groupwidth / (2*numbars);
    end
    handles.errors(b) = errorbar(x', meanSideBySide(:,b), stdSideBySide(:,b), 'k', 'linestyle', 'none', 'linewidth', lw);
end

colormap 'summer'; % axis tight; % should be inside plotopt

if isfield(plotopt, 'labels')
    set(gca,'XTickLabel', plotopt.labels);
    set(gca,'XLim', [0 numel(plotopt.labels)+2])
    set(gca,'XTick', (1:numel(plotopt.labels)))
    rotateticklabel(gca, 45);
end
if isfield(plotopt, 'title')
    title(plotopt.title, 'FontSize', fn);
end
if isfield(plotopt, 'ylabel')
    ylabel(plotopt.ylabel, 'FontSize', fn);
end
if isfield(plotopt, 'xlabel')
    xlabel(plotopt.xlabel, 'FontSize', fn);
end
if isfield(plotopt, 'legend')
    legend(plotopt.legend, 'Location', 'Best');
end
