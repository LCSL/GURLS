function [] = summary_plot(filestr, fields, nRuns, plotopt, fileroot)
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
    for indFile = 1:numel(filestr)
        meanSideBySide(:, indFile) = res{indFile}{f}.perclassMean';
        stdSideBySide(:, indFile) = res{indFile}{f}.perclassStd';
    end
    
    if isfield(plotopt, 'titles')
        op.title = plotopt.titles{f};
    end
    if isfield(plotopt, 'ylabels')
        op.ylabel = plotopt.ylabels{f};
    end
    if isfield(plotopt, 'xlabels')
        op.xlabel = plotopt.xlabels{f};
    end
    if isfield(plotopt, 'class_names')
        op.labels = plotopt.class_names;
    end
    if isfield(plotopt, 'legend')
        op.legend = plotopt.legend;
    else
        op.legend = filestr;
    end
    
    pretty_plot(meanSideBySide, stdSideBySide, op);

end

