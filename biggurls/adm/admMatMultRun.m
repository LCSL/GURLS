function [out] = admMatMultRun(jobFileName)
		jf = load(jobFileName);
		X = bigarray.Obj(jf.jobStruct.XPath);
		Y = bigarray.Obj(jf.jobStruct.YPath);
		fLock = admSetup(jf.jobStruct.stateFileName);
		
		blockID = NO_WORK;
		
		while (blockID ~= JOB_DONE)
				blockID = getWork(jf.jobStruct.stateFileName, fLock);
				if ~exist(jf.jobStruct.stateFileName,'file')
					fprintf('Unable to read "%s"\n', jf.jobStruct.stateFileName);
					blockID = JOB_DONE;
				end	
				switch blockID
					case JOB_DONE
						fprintf('Job Completed!\r');
					case NO_WORK
						fprintf('Node out of work\r');
						pause(10);
					case FINALIZE_JOB
						fprintf('Finalizing Job\r');
						finalizeJob(jf.jobStruct.stateFileName, X.NumBlocks, jobFileName);
						reportWork(jf.jobStruct.stateFileName, fLock, blockID);
					otherwise
						fprintf('Taking care of job: %d\r', blockID);
						d = multBlock(X,Y,double(blockID));
						save([jf.jobStruct.stateFileName '_' num2str(blockID)], 'd');
						reportWork(jf.jobStruct.stateFileName, fLock, blockID);
				end
		end
		
		admDismiss(fLock);
		if exist(jf.jobStruct.stateFileName,'file')
			delete(jf.jobStruct.stateFileName);
		end	
end		
function d = multBlock(X,Y,blockID)
		bX = X.ReadBlock(blockID);
		bY = Y.ReadBlock(blockID);

		d = bX*bY';
end

function [] = finalizeJob(stateFileName, nB, jobFileName)
		t = load([stateFileName '_' num2str(1)]);
		data = zeros(size(t.d));
		for i = 1:nB
			bName = [stateFileName '_' num2str(i)];
			t = load(bName);
			delete([bName '.mat']);
			data = data + t.d;
		end
		save(jobFileName, 'data');
end	



