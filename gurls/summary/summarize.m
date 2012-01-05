function [res] = summarize(filestr, fields, nRuns)
% Give filenames
% field to summarize.
% it is user's responsibility to make sure the said feild exists in all the files, in the right format,
% and can be averaged.
% The average is computed over the files for the field.



for run = 1:nRuns,
    try
    file = ([filestr '_' num2str(run)]);
    load(file);
    catch % This is to understand to neva/jim's convention, which by the way must be stopped!
    file = sprintf('%s_%02s', filestr, num2str(run));
    load(file);
    end
    fprintf('Run: ..%d ', run);
    for f = 1:numel(fields)
        res{f}.name = fields{f};
        s = regexp(fields{f},'\.','split');
        % Is copying to tmp necessary? This line may be deleted!
        tmp = opt;
        % This if condition is added for backward compatiability; 
        % with the previous version where perf was stored as an array
        % of structures, now it is stored as stucture of arrays. 
        if numel(tmp.(s{1})) > 1
            for t = 1:numel(tmp.(s{1}))
                res{f}.val(run,t) = getfield(opt.(s{1})(t), s{2});
            end
            
        else
            
            for i = 1:numel(s)
                tmp = getfield(tmp,s{i});
            end
            res{f}.val(run,:) = tmp;
        end
    end
end
fprintf('done.\n');
for f = 1:numel(fields)
    res{f}.perclassMean = mean(res{f}.val,1);
    res{f}.perclassStd  = std (res{f}.val,0,1);
    
    overall = median(res{f}.val'); % Should this be an option? @mean/@median
    res{f}.overallMean = mean(overall);
    res{f}.overallStd  = std(overall);
end
