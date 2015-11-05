function res = summarize(filestr, fields, nRuns)
% Compute average performance metric over all .mat files (gurl outputs)
% for a given field.
%
% NOTE: it is user's responsibility to make sure that 'fields' exists in
% all files/trials descibed by filestr and nRuns, in the right format, and
% can be averaged.
%
% INPUTS:
% - filestr: basename (str) of file (.mat) to load
% - fields: field in loaded structure to summarize.
% - nRuns: number of runs, i.e. filenames of the form [sprintf(filestr_%d, nRuns)]
%
% OUTPUT:
% - res: a struct with summary fields, e.g. perclassMean, perclassStd, overallMean, overallStd:

nFields = numel(fields);

for r = 1:nRuns,
    
    fprintf('Run: %d \n', r);
    try
        file = ([filestr '_' num2str(r)]);
        load(file);
    catch % Take care of older naming convention which is now obsolete!
        file = sprintf('%s_%02s', filestr, num2str(r));
        load(file);
    end
    
    for f = 1:nFields
        res{f}.name = fields{f};
        s = regexp(fields{f}, '\.', 'split');
        
        % This if condition is added for backward compatibility;
        % with the previous version where perf was stored as an array
        % of structures, now it is stored as stucture of arrays.
        
        opt_tmp = opt; % Is copying to tmp necessary?
        if numel(opt_tmp.(s{1})) > 1
            for t = 1:numel(opt_tmp.(s{1}))
                res{f}.val(r,t) = opt.(s{1})(t).(s{2});
            end
        else
            % weird recursion to get the perf value
            for i = 1:numel(s)
                opt_tmp = opt_tmp.(s{i});
            end
            res{f}.val(r,:) = opt_tmp;
        end
    end
end
fprintf('done.\n');

for f = 1:nFields
    res{f}.perclassMean = mean(res{f}.val, 1);
    res{f}.perclassStd = std(res{f}.val, 0, 1);
    
    overall = median(res{f}.val, 2)'; % Should this be an option? @mean/@median
    res{f}.overallMean = mean(overall);
    res{f}.overallStd = std(overall);
end
