#include "mex.h"
#include <unistd.h>
#include "../include/lockLib.h"
#include "../include/retCode.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		if(nrhs != 3) mexErrMsgTxt("reportWork: incorrect number of arguments");
		if(nlhs != 0) mexErrMsgTxt("reportWork: incorrect number of outputs");

		char stateFileName[512];
		int fLock;
		int finishedBlockID;

		mxGetString(prhs[0], stateFileName, 512);
		fLock = (int) mxGetScalar(prhs[1]);
		finishedBlockID = (int) mxGetScalar(prhs[2]);

		while(reportWork(stateFileName, fLock, finishedBlockID) != REPORT_SUCCESS) sleep(1);
}
