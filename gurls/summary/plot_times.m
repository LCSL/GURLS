function [] = plot_times(filestr, nRuns, res_root, legend_str)

if ~exist('res_root', 'var')
    res_root = '.';
end

u_fieldlist = {};

for indFile = 1:numel(filestr)
    
    finalStruct{indFile} = struct;
    
    for r = 1:nRuns(indFile)
        
        fprintf('Run: %d \n', r);
        try
            file = ([filestr{indFile} '_' num2str(r)]);
            load(fullfile(res_root, file));
        catch % Take care of older naming convention which is now obsolete!
            file = sprintf('%s_%02s', filestr{indFile}, num2str(r));
            load(fullfile(res_root, file));
        end
        
        numJobs = numel(opt.time);
        for job = 1:numJobs
            fieldlist = fieldnames(opt.time{job});
            for field = 1:numel(fieldlist)
                if ~isfield(finalStruct{indFile},fieldlist{field})
                    finalStruct{indFile}.(fieldlist{field}) = [];
                end
                curr = finalStruct{indFile}.(fieldlist{field});
                v = opt.time{job}.(fieldlist{field});
                finalStruct{indFile}.(fieldlist{field}) = [curr v];
            end
        end
    end
    finalStruct{indFile}.total = sum(cell2mat(struct2cell(finalStruct{indFile})));
end
fprintf('done.\n');

for indFile = 1:numel(filestr)
    fieldlist = fieldnames(finalStruct{indFile});
    meanStruct{indFile} = struct;
    stdStruct{indFile} = struct;
    for field = 1:numel(fieldlist)
        s = std(finalStruct{indFile}.(fieldlist{field}));
        m = mean(finalStruct{indFile}.(fieldlist{field}));
        stdStruct{indFile}.(fieldlist{field}) = s;
        meanStruct{indFile}.(fieldlist{field}) = m;
    end
    u_fieldlist = union(u_fieldlist, fieldlist);
end
for field = 1:numel(u_fieldlist)
    for indFile = 1:numel(filestr)
        if ~isfield(meanStruct{indFile},u_fieldlist{field});
            stdStruct{indFile}.(u_fieldlist{field}) = 0;
            meanStruct{indFile}.(u_fieldlist{field}) = 0;
        end
        meanSideBySide(field,indFile) = meanStruct{indFile}.(u_fieldlist{field});
        stdSideBySide(field,indFile) = stdStruct{indFile}.(u_fieldlist{field});
    end
end

%%%% PLOT %%%%
plotopt.title = 'Time performance comparison';
plotopt.ylabel = '[s]';
plotopt.labels = u_fieldlist;
if nargin < 4
    plotopt.legend = filestr;
else
    plotopt.legend = legend_str;
end
pretty_plot(meanSideBySide, stdSideBySide, plotopt);

