#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include <math.h>

#ifndef max
	#define max(x,y)      ((x) < (y) ? (y) : (x))
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *K, *u, *v, *wj, *activeset;
  double one = 1.0, zero = 0.0;
  int oneint = 1; 
  int N,M,L;
  int j,ind;
  int *Kdims;
  bool b_active;
  char *chn = "N";
  

  K = mxGetPr(prhs[0]);
  u = mxGetPr(prhs[1]);
 
  Kdims = (int*)mxGetDimensions(prhs[0]);
  N = Kdims[0];
  /*//N2 = Kdims[0];*/
  M = Kdims[2];
  
  if(nrhs >= 3){
      b_active = 1;
      L = (int)mxGetN(prhs[2]);
      activeset = mxGetPr(prhs[2]);
  }else{
      b_active = 0;
  };
  

  /*//plhs[0] = (mxArray *)prhs[2];
  //plhs[1] = (mxArray *)prhs[3]; */

    if(b_active){ 
        plhs[0] = mxCreateDoubleMatrix(N,L,mxREAL);
    }else{
        plhs[0] = mxCreateDoubleMatrix(N,M,mxREAL);
    };
  v = mxGetPr(plhs[0]);

  if(nlhs>=2){
      if(b_active){ 
          plhs[1] = mxCreateDoubleMatrix(L,1,mxREAL);
      }else{
          plhs[1] = mxCreateDoubleMatrix(M,1,mxREAL);
      };
      wj = mxGetPr(plhs[1]);
  };
  
  if(b_active){
      #pragma omp parallel for
      for(j=0;j<L;j++){
          ind = (int)activeset[j]-1;
          dgemv(chn, &N, &N, &one, &K[ind*N*N], &N, &u[ind*N], &oneint, &zero, &v[j*N], &oneint);
          if(nlhs>=2){
              wj[j] = sqrt(ddot(&N,&u[ind*N], &oneint, &v[j*N], &oneint));              
          };
          
      };
  }else{
      #pragma omp parallel for
      for(j=0;j<M;j++){
          dgemv(chn, &N, &N, &one, &K[j*N*N], &N, &u[j*N], &oneint, &zero, &v[j*N], &oneint);
          if(nlhs>=2){
              wj[j] = sqrt(ddot(&N,&u[j*N], &oneint, &v[j*N], &oneint));
          };

      };
  };   
  return;
  
};

