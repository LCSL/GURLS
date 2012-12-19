#include "mex.h"
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

#include "../include/retCode.h"
#include "../include/lockLib.h"

int getWork(const char *stateFileName, int fLock)
{
		FILE *fState;

		if(acquireLock(fLock) == LOCK_FAIL) return NO_WORK;

		if((fState = fopen(stateFileName, "r+")) == NULL)
		{
				mexPrintf("Unable to read stateFile\n");
				fState = NULL;
				return JOB_DONE;
		}

		int blockID;
		int state;
		int dummy;
		int curPos;

		while (!feof(fState))
		{
				dummy = fscanf(fState,"%d\t%d\n",&blockID, &state);
				if(state == 0) break;
		}

		/*
		 *  Normal case: we found a new job to do.
		 *	We just want to check whether it is a normal 
		 *	job or a finalizing job.
		 */

		//mexPrintf("blockID:\t%d\nstate:\t\t%d\n",blockID,state);

		switch(state)
		{
				case 0:
						if(blockID > 0)
						{
								/*
								 *	Normal case:
								 *	- regular job found
								 *	- tell the world we'll take care of it
								 */
								curPos = ftell(fState);
								fseek(fState, curPos - 2, SEEK_SET);
								fputs("1\n",fState);
						}
						else
						{
								/*
								 * 	Finalize task:
								 * 	- the first available job is the last one (finalization)
								 *	- make sure all other jobs have been completed and not jus booked
								 *	- tell the workd we'll take care of this task
								 */

								curPos = ftell(fState);
								fseek(fState, 0, SEEK_SET);
								while(!feof(fState))
								{
										dummy = fscanf(fState,"%d\t%d\n",&blockID, &state);
										if(state != 2 && blockID != -1)
										{
												blockID = NO_WORK;
												break;
										}
								}
								if(blockID == FINALIZE_JOB)
								{
										fseek(fState, curPos - 2, SEEK_SET);
										fputs("1\n",fState);
								}
						}
						break;
				case 2:
						blockID = JOB_DONE;
						break;
				default:
						blockID = NO_WORK;
						break;
		}

		fclose(fState);
		while(releaseLock(fLock) != LOCK_SUCCESS) sleep(1);
		return blockID;
}

int reportWork(const char *stateFileName, int fLock, int finishedBlockID)
{
		FILE *fState;

		if(acquireLock(fLock) == LOCK_FAIL) return REPORT_FAIL;

		if((fState = fopen(stateFileName, "r+")) == NULL)
		{
				mexPrintf("Unable to read stateFile\n");
				fState = NULL;
				return REPORT_FAIL;
		}

		int blockID;
		int state;
		int dummy;


		while (!feof(fState))
		{
				dummy = fscanf(fState,"%d\t%d\n",&blockID, &state);
				if(blockID == finishedBlockID) break;
		}


		int curPos = ftell(fState);
		fseek(fState, curPos - 2, SEEK_SET);
		fputs("2\n",fState);
		fclose(fState);

		while(releaseLock(fLock) != LOCK_SUCCESS) sleep(1);
		return REPORT_SUCCESS;
}
