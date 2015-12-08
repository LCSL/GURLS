#include "mex.h"
#include <math.h>

/*
 * This function computes
 * f = - 0.5 * \sum_{i,j}\alpha_i\alpha_jy_iy_j exp(-\sum_k theta(k)*d_k(x_i,x_j))
 * df(l) = 0.5 * \sum_{i,j}\alpha_i\alpha_jy_iy_j d_l(x_i,x_j) *exp(-\sum_k theta(k)*d_k(x_i,x_j))
 * H(k1,k2) = -0.5 * \sum_{i,j}\alpha_i\alpha_jy_iy_j d_k2(x_i,x_j) d_k1(x_i,x_j) *exp(-\sum_k theta(k)*d_k(x_i,x_j))
 *
 * input 
 *		D : [n,n,nFactors]
 *		fac: [nFactors,1]
 *		alpha: [n,1]
 *
 * output
 *		f = [1,1];
 *		df = [nFactors,1];
 *		H = [nFactors,nFactors]
 *
 * important!!!! This function assumes d_ii = 0 and d_ij=d_ji
 */    
 void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if ((nrhs!=3) & (nrhs!=4))
		mexErrMsgTxt("Wrong number of input arguments");

	if ((nlhs==0) |(nlhs>3))
		mexErrMsgTxt("Wrong number of output arguments min 1, max 3");

	if (mxGetNumberOfDimensions(prhs[1])!=2)
		mexErrMsgTxt("Second input argument must be 2D");

	const unsigned int nDims = mxGetNumberOfDimensions(prhs[0]);
	const mwSize *dims = mxGetDimensions(prhs[0]);

	const unsigned int nFactors = (nDims!=3)?1:dims[2];

	assert(nFactors>0);

	const double *Ds = mxGetPr(prhs[0]);
	const double *factors = mxGetPr(prhs[1]);
	const double *alpha = mxGetPr(prhs[2]);

	const bool smoothing = (nrhs==4)?true:false;

	const unsigned int n = mxGetM(prhs[2]);

	if ((n != dims[1]) || (n!=dims[0]))
		mexErrMsgTxt("First argument must be symmetric in first 2 dimensions");
	
	const unsigned int nn = n*n;
	
	if (mxGetM(prhs[1]) != nFactors)
		mexErrMsgTxt("Second argument must be of size [Nfactors,1]");
	
	if (mxGetN(prhs[1]) != 1)
		mexErrMsgTxt("Second argument must be of size [Nfactors,1]");

	if (mxGetN(prhs[2]) != 1)
		mexErrMsgTxt("Third argument must be of size [Npoints,1]");

	// output arguments
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
	double *f = mxGetPr(plhs[0]); // objective value

	double *df = NULL; // derivative
	if (nlhs>1) {
		plhs[1]=mxCreateDoubleMatrix(nFactors,1,mxREAL); 
		df = mxGetPr(plhs[1]);
	}
	
	double *H = NULL; // Hessian
	if (nlhs>2) {
		plhs[2]=mxCreateDoubleMatrix(nFactors,nFactors,mxREAL); 
		H = mxGetPr(plhs[2]);
	}

	double scale = 1.0;
	if (!smoothing)
	{
		for (unsigned int i=0; i<n-1 ; ++i)
		{
			for (unsigned int j=i+1,jn=(i+1)*n; j<n ; ++j,jn+=n)
			{

				double d = 0.0;
				const unsigned int ijn = i+jn;
			
				for (unsigned int k=0 ; k<nFactors; ++k)
					d += factors[k] * Ds[ijn+k*nn];
				const double w = alpha[i]*alpha[j]*exp(-d);

				*f -= w;

				if (df)
					for (unsigned int k=0 ; k<nFactors ; ++k)
					{
						const double tmp = w*Ds[ijn+k*nn];
						df[k] += tmp;
						if (H)
							for (unsigned int l=0 ; l<=k ; ++l)
								H[k+l*nFactors] -= tmp*Ds[ijn+l*nn];
					}
			}
		}
	}
	else
	{
		double *sigma = mxGetPr(prhs[3]);
		if (*sigma<0)
			mexErrMsgTxt("Smoothing variance must be positive");

		scale = pow(sigma[0] * (3.141593),(-((double)nFactors+1)/2));
		for (unsigned int i=0; i<n-1 ; ++i)
		{
			for (unsigned int j=i+1,jn=(i+1)*n; j<n ; ++j,jn+=n)
			{

				double d = 0.0;
				double s = 1.0;
				for (unsigned int k=0 ; k<nFactors; ++k)
				{
					d += factors[k] * (Ds[i+jn+k*nn]/(1+Ds[i+jn+k*nn]*sigma[0]));
					s /= sqrt(1+sigma[0]*Ds[i+jn+k*nn]);
				}

				const double w = alpha[i]*alpha[j]*s*exp(-d);

				*f -= w;

				if (df)
					for (unsigned int k=0 ; k<nFactors ; ++k)
					{
						const double tmp = w*Ds[i+jn+k*nn];
						df[k] += tmp;
						if (H)
							for (unsigned int l=0 ; l<=k ; ++l)
								H[k+l*nFactors] -= tmp*Ds[i+jn+l*nn];
					}
			}
		}
		
	}

	
	double alpha2 = 0.0;
	for (unsigned int i=0 ; i<n ; ++i)
		alpha2 += alpha[i]*alpha[i];

	*f -= 0.5*alpha2;


	// need to correct the range of the objective function

/*
  if (smoothing)
  {
  *f *= (scale);
  if (df)
  for (unsigned int k=0 ; k<nFactors ; ++k)
  {
  df[k] *= scale;
  if (H)
  for (unsigned int l=0 ; l<=k ; ++l)
  H[k+l*nFactors] *= scale;
  }
  }
*/

}


// If Matlab matrix is 3D A, then A(r,c,d)
// pSize = mxGetDimensions(prhs[0]);
// pA = mxGetPr(prhs[0]);
// A_rcd = pA[r-1 + (c-1)*pSize[0] + (d-1)*pSize[0]*pSize[1]];
