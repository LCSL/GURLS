#include "mex.h"
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

#include "../include/lockLib.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		if(nrhs != 1) mexErrMsgTxt("getWork: incorrect number of arguments");
		if(nlhs != 0) mexErrMsgTxt("getWork: incorrect number of outputs");

		int fLock;

		fLock = (int) mxGetScalar(prhs[0]);
		close(fLock);
}


