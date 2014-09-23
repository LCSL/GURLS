 /*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
  * All rights reserved.
  *
  * authors:  M. Santoro
  * email:   msantoro@mit.edu
  * website: http://cbcl.mit.edu/IIT@MIT/IIT@MIT.html
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions
  * are met:
  *
  *     * Redistributions of source code must retain the above
  *       copyright notice, this list of conditions and the following
  *       disclaimer.
  *     * Redistributions in binary form must reproduce the above
  *       copyright notice, this list of conditions and the following
  *       disclaimer in the documentation and/or other materials
  *       provided with the distribution.
  *     * Neither the name(s) of the copyright holders nor the names
  *       of its contributors or of the Massacusetts Institute of
  *       Technology or of the Italian Institute of Technology may be
  *       used to endorse or promote products derived from this software
  *       without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGE.
  */


#ifndef _GURLS_RLSPRIMALRECUPDATECHOLESKY_H_
#define _GURLS_RLSPRIMALRECUPDATECHOLESKY_H_

#include "gurls++/optimization.h"

#include "gurls++/optmatrix.h"
#include "gurls++/optfunction.h"

#include "gurls++/utils.h"
#include "gurls++/cholupdateutils.h"    // Include functions for Cholesky update via BLAS calls
#include "gmat2d.h"


namespace gurls
{

/**
 * \ingroup Optimization
 * \brief RLSPrimalRecUpdateCholesky is the sub-class of Optimizer that updates the upper Cholesky factor of the covariance matrix
 */
template <typename T>
class RLSPrimalRecUpdateCholesky: public Optimizer<T>
{
public:
    /**
     * Computes a classifier for the primal formulation of RLS, using a
     * recursive Cholesky update, starting from an initial estimator found in opt.optimizer.
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following fields that need to be set through previous gurls++ tasks:
     *  - optimizer.W (settable with the class RLSPrimalRecInitCholesky)
     *  - optimizer.R (settable with the class RLSPrimalRecInitCholesky)
     *  - optimizer.b (settable with the class RLSPrimalRecInitCholesky)
     *
     * \return adds to opt the field optimizer which is a list containing the following fields:
     *  - W = updated matrix of coefficient vectors of rls estimator for each class
     *  - R = Updated Cholesky factor of the covariance matrix
     *  - b = Updated Xty matrix
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* RLSPrimalRecUpdateCholesky<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{   
    // Get significant sizes
    //	[n,d] = size(X);

    const unsigned long n = X.rows();   // Number of updates to be performed
    const unsigned long d = X.cols();
    const unsigned long t = Y.cols();

    // Xty = X'*y;
    T* Xty;
    Xty = new T[d*t];
    dot(X.getData(), Y.getData(), Xty, n, d, n, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
    
    // Retrieve previous W, R, b
    
    //  W = opt.rls.W;
    const gMat2D<T>& prev_W = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.W");
    gMat2D<T>* W = new gMat2D<T>(prev_W);
    
    //  R = opt.rls.R;
    const gMat2D<T>& prev_R = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.R");
    gMat2D<T>* R = new gMat2D<T>(prev_R);
    
    //  b = opt.rls.b;
    const gMat2D<T>& prev_b = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.b");
    gMat2D<T>* b = new gMat2D<T>(prev_b);
    
    // Update R
    for(unsigned long i=0; i<n; i++)
    {
        // Get pointer to the i-th input
//         gVec<T> v( X[i] );
//         double* input =  v.getData();
        T* input = (T*)X.getData() + i;
        
        // Rank-1 Cholesky update of R
        cholupdate(*R, *input, 1);

        // Update b = b + X[i]ty
        // TODO: Maybe the update of b can be done in 1-shot outside the for loop.
        gMat2D<T> transposedinput(d, 1);
        X.transpose(transposedinput);
        gMat2D<T> res( d , t );          // Temporary variable for result storage
        
        dot( transposedinput , Y , res );
        *b = (*b) + res;
    }

    // Update W
    
    copy(W->getData(), b->getData(), W->getSize());
    
    // Forward substitution
    mldivide_squared(R->getData(), W->getData(), d, d, W->rows(), W->cols(), CblasTrans);    //(R'\Xty) 
    
    // Backward substitution
    mldivide_squared(R->getData(), W->getData(), d, d, W->rows(), W->cols(), CblasNoTrans);  //R\(R'\Xty)     

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

    // cfr.R = R;
    optimizer->addOpt("R", new OptMatrix<gMat2D<T> >(*R));
    
    //  cfr.W = W;
    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    //  cfr.b = b;
    optimizer->addOpt("b", new OptMatrix<gMat2D<T> >(*b));
    
    delete[] Xty;
    
    return optimizer;
}


}
#endif // _GURLS_RLSPRIMALRECUPDATE_H_
