 /*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
  * All rights reserved.
  *
  * author:  Raffaello Camoriano
  * email:   raffaello.camoriano@iit.it
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


#ifndef _GURLS_RLSPRIMALRECINITCHOLESKY_H_
#define _GURLS_RLSPRIMALRECINITCHOLESKY_H_

#include "gurls++/optimization.h"

#include "gurls++/optmatrix.h"
#include "gurls++/optfunction.h"

#include "gurls++/utils.h"

namespace gurls
{

/**
 * \ingroup Optimization
 * \brief RLSPrimalRecInitCholesky is the sub-class of Optimizer that implements RLS with the primal formulation and computes the upper Cholesky factor of the regularized covariance matrix.
 */
template <typename T>
class RLSPrimalRecInitCholesky: public Optimizer<T>
{

public:
    /**
     * Computes a classifier for the primal formulation of RLS, storing the Cholesky decomposition of the 
     * covariance matrix.
     * The regularization parameter is set to the one found in the field paramsel of opt.
     * The variables necessary for further recursive update (R, W, b) are stored in the output structure
     * In case of multiclass problems, the regularizers need to be combined with the function specified inthe field singlelambda of opt
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - singlelambda (default)
     *  - paramsel (settable with the class ParamSelection and its subclasses)
     *
     * \return adds to opt the field optimizer which is a list containing the following fields:
     *  - W = matrix of coefficient vectors of rls estimator for each class
     *  - R = Upper Cholesky factor of the regularized covariance matrix (XtX + n*lambda*I)
     *  - b = Xty vector
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* RLSPrimalRecInitCholesky<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //lambda = opt.singlelambda(opt.paramsel.lambdas);
    const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
    T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());
    
    // NOTE: An error should be thrown if paramsel.lambdas is not defined

    //[n,d] = size(X);
    const unsigned long n = opt.hasOpt("nTot")? static_cast<unsigned long>(opt.getOptAsNumber("nTot")) : X.rows();
    unsigned long d;
    unsigned long t;

    //XtX = X'*X;
    T* XtX;
    if(!opt.hasOpt("kernel.XtX"))
    {
        d = X.cols();
        XtX = new T[d*d];
        dot(X.getData(), X.getData(), XtX, n, d, n, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
    }
    else
    {
        const gMat2D<T>& XtX_mat = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.XtX");
        d = XtX_mat.cols();
        XtX = new T[d*d];
        copy(XtX, XtX_mat.getData(), d*d);
    }

    //Xty = X'*y;
    T* Xty;
    if(!opt.hasOpt("kernel.Xty"))
    {
        t = Y.cols();
        Xty = new T[d*t];
        dot(X.getData(), Y.getData(), Xty, n, d, n, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
    }
    else
    {
        const gMat2D<T>& Xty_mat = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.Xty");
        t = Xty_mat.cols();
        Xty = new T[d*t];
        copy(Xty, Xty_mat.getData(), d*t);
    }

    // Compute R starting from XtX
    // R = chol(XtX + (n*lambda)*eye(d));

    T coeff = n*lambda;
    axpy(d, (T)1.0, &coeff, 0, XtX, d+1);	// Regularize XtX

    // Cholesky factorization

    //R = chol(XtX);
    T* R = new T[d*d];
    cholesky(XtX, d, d, R);     // Computes the upper Cholesky factor by default
    gMat2D<T> *R_mat = new gMat2D<T>(R, d, d, 1);
    
    // Initialize weights vector
    // W = R\(R'\Xty);
    gMat2D<T> *W = new gMat2D<T>(d, t);

    copy(W->getData(), Xty, W->getSize());
    
    // Forward substitution
    mldivide_squared(R, W->getData(), d, d, W->rows(), W->cols(), CblasTrans);    //(R'\Xty)
    
    // Backward substitution
    mldivide_squared(R, W->getData(), d, d, W->rows(), W->cols(), CblasNoTrans);  //R\(R'\Xty)
    
    // Initialize b
    gMat2D<T> *b = new gMat2D<T>(Xty, d, t, 1);
    
    delete[] XtX;
    
    //Save results in OPT
    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");
    
    std::cout << "W:" << std::endl;
    std::cout << *W << std::endl;
    std::cout << "b:" << std::endl;
    std::cout << *b << std::endl;
    std::cout << "R:" << std::endl;
    std::cout << *R_mat << std::endl;
    
    //  cfr.W = W;
    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    //  cfr.b = b;
    optimizer->addOpt("b", new OptMatrix<gMat2D<T> >(*b));

    //  cfr.R = R;
    optimizer->addOpt("R", new OptMatrix<gMat2D<T> >(*R_mat));

    
    return optimizer;
}


}
#endif // _GURLS_RLSPRIMALRECINITCHOLESKY_H_
