#include "mex.h"
#include "../include/lockLib.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		if(nrhs != 2) mexErrMsgTxt("getWork: incorrect number of arguments");
		if(nlhs != 1) mexErrMsgTxt("getWork: incorrect number of outputs");

		char stateFileName[512];
		int fLock;
		int blockID;

		mxGetString(prhs[0], stateFileName, 512);
		fLock = (int) mxGetScalar(prhs[1]);

		blockID = getWork((const char *) stateFileName, fLock);
		plhs[0] = mxCreateNumericMatrix(1,1,mxINT64_CLASS, mxREAL);
		*(int64_T*) mxGetData(plhs[0]) = (int64_T) blockID;
}

