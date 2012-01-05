function [] = plot_times(filestr,nRuns, legend_str)

u_fieldlist = {};		

for experiment = 1:numel(filestr)
	finalStruct{experiment} = struct;	
	for run = 1:nRuns{experiment}
        fprintf('Run: ..%d ', run);
        try 
		file = ([filestr{experiment} '_' num2str(run)]);
        load(file);
        catch % This is to confirm to Neva/Jim's filenaming;
        file = sprintf('%s_%02d',filestr{experiment},run);
		load(file);
        end
		numJobs = numel(opt.time);
		for job = 1:numJobs	
			fieldlist = fieldnames(opt.time{job});
			for field = 1:numel(fieldlist)
				if ~isfield(finalStruct{experiment},fieldlist{field})
					finalStruct{experiment} = setfield(finalStruct{experiment},fieldlist{field},[]);
				end
				curr = getfield(finalStruct{experiment},fieldlist{field});
				v = getfield(opt.time{job},fieldlist{field});
				finalStruct{experiment} = setfield(finalStruct{experiment},fieldlist{field},[curr v]);	
			end
		end	
	end
	finalStruct{experiment}.total = sum(cell2mat(struct2cell(finalStruct{experiment})));
end	
fprintf('done.\n');

for experiment = 1:numel(filestr)
	fieldlist = fieldnames(finalStruct{experiment});	
	meanStruct{experiment} = struct;
	stdStruct{experiment} = struct;
	for field = 1:numel(fieldlist)	
		s = std(getfield(finalStruct{experiment},fieldlist{field}));
		m = mean(getfield(finalStruct{experiment},fieldlist{field}));
		stdStruct{experiment} = setfield(stdStruct{experiment},fieldlist{field},s);
		meanStruct{experiment} = setfield(meanStruct{experiment},fieldlist{field},m);
	end
	u_fieldlist = union(u_fieldlist, fieldlist);
end	
for field = 1:numel(u_fieldlist)
	for experiment = 1:numel(filestr)
		if ~isfield(meanStruct{experiment},u_fieldlist{field});
			stdStruct{experiment} = setfield(stdStruct{experiment},u_fieldlist{field},0);
			meanStruct{experiment} = setfield(meanStruct{experiment},u_fieldlist{field},0);
		end	
		meanSideBySide(field,experiment) = getfield(meanStruct{experiment},u_fieldlist{field});
		stdSideBySide(field,experiment) = getfield(stdStruct{experiment},u_fieldlist{field});
	end
end		

%%%% PLOT %%%%
plotopt.title = 'Time performance comparison';
plotopt.ylabel = '[s]';
plotopt.labels = u_fieldlist;
if nargin < 3
	plotopt.legend = filestr;
else
    plotopt.legend = legend_str;
end
pretty_plot(meanSideBySide,stdSideBySide,plotopt);

