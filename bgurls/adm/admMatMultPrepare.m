function [] = admMatMultPrepare(X,Y, jobFileName)
		jobStruct = struct;
		jobStruct.XPath = X.Path;
		jobStruct.YPath = Y.Path;
		jobStruct.stateFileName = [pwd() '/' datestr(now,30)];

		f = fopen(jobStruct.stateFileName, 'w');
		for s = 1:X.NumBlocks
			fprintf(f,'%d\t%d\n',s,0);
		end
		fprintf(f,'%d\t%d\n',-1,0);
		fclose(f);
		save(jobFileName,'jobStruct');
		pause(5);
