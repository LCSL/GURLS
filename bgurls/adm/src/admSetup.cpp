#include "mex.h"
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

#include "../include/lockLib.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		if(nrhs != 1) mexErrMsgTxt("getWork: incorrect number of arguments");
		if(nlhs != 1) mexErrMsgTxt("getWork: incorrect number of outputs");

		char stateFileName[512];
		int fLock;
		mxGetString(prhs[0], stateFileName, 512);
		fLock = open(stateFileName, O_RDWR);
		if(fLock == -1)
		{
				mexErrMsgTxt("Could not set lock on file");
		}
		
		plhs[0] = mxCreateNumericMatrix(1,1,mxINT64_CLASS, mxREAL);
		*(int64_T*) mxGetData(plhs[0]) = (int64_T) fLock;
}


