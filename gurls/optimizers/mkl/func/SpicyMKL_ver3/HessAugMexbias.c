#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include <math.h>
#include <string.h>

/*#include "routines.h"*/

#ifdef _OPENMP
#include <omp.h>
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *K, *activeset, *v, *wj,  *cgamma, *HessAug, *alpha1, *alpha2;
    int ind;
    double one = 1.0, zero = 0.0, alpha, C, cgammab, work;
    int oneint = 1; 
    int N,M,len, NN;
    int j,i,k,l,m;
    int *Kdims;
    char *chn = "N";

    K = mxGetPr(prhs[0]);
    wj = mxGetPr(prhs[1]);
    activeset = mxGetPr(prhs[2]);
    v = mxGetPr(prhs[3]);
    C = mxGetScalar(prhs[4]);
    cgamma = mxGetPr(prhs[5]);
    cgammab = mxGetScalar(prhs[6]);
    alpha1 = mxGetPr(prhs[8]);
    alpha2 = mxGetPr(prhs[9]);
    

    Kdims = (int*)mxGetDimensions(prhs[0]);
    N = Kdims[0];
    M = Kdims[2];
    len = mxGetN(prhs[2]);
    NN = N*N;

    plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
    HessAug = mxGetPr(plhs[0]);
    memcpy(HessAug,mxGetPr(prhs[7]),N*N*sizeof(HessAug[0]));
    
  
#ifdef _OPENMP
      if(NN > 10000){      
          omp_set_num_threads(4);
      }else{
          omp_set_num_threads(2);
      }

    for(j=0;j<len;j++){
        ind = (int)activeset[j] -1;
        
        m = ind*NN;
        #pragma omp parallel for private(i) shared(HessAug) 
        for(i=0;i<NN;i++){
            HessAug[i] += alpha1[ind]*K[m + i]; 
        };
        m = ind*N;
        #pragma omp parallel for private(i,l,k) shared(HessAug) 
        for(i=0;i<N;i++){
          l = i*N;
          work = alpha2[ind]*v[m + i];
          for(k=0;k<N;k++){
              HessAug[l + k] += v[m + k]*work;
          };
        };
    };
    #pragma omp flush
#else 
    for(j=0;j<len;j++){
        ind = (int)activeset[j] -1;        
        daxpy(&NN, &alpha1[ind], &K[NN*ind], &oneint, HessAug, &oneint);
        dger(&N, &N, &alpha2[ind], &v[ind*N], &oneint, &v[ind*N], &oneint, HessAug, &N );
    };
#endif
    
    #pragma omp parallel for private(j)
    for(j=0;j<NN;j++){
        HessAug[j] += cgammab;   
    };
    return;
};
