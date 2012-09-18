 /*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011, IIT@MIT Lab
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


#ifndef _GURLS_RLSPRIMALR_H_
#define _GURLS_RLSPRIMALR_H_

#include "optimization.h"
#include "gmath.h"
#include "utils.h"

namespace gurls {

/**
 * \ingroup Optimization
 * \brief RLSPrimalr is the sub-class of Optimizer that implements RLS with the primal formulation, using a randomized version of SVD
 */
template <typename T>
class RLSPrimalr: public Optimizer<T>{

public:
    /**
     * Computes a classifier for the primal formulation of RLS, using a randomized version of Singular value decomposition.
     * The regularization parameter is set to the one found in the field paramsel of opt.
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
     *  - C = empty matrix
     *  - X = empty matrix
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* RLSPrimalr<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)
{
    //	lambda = opt.singlelambda(opt.paramsel.lambdas);
    const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
    T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());


//    std::cout << "Solving primal RLS using Randomized SVD..." << std::endl;

    //	[n,d] = size(X);

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();

    const unsigned long Yn = Y.rows();
    const unsigned long Yd = Y.cols();

    //	===================================== Primal K

//	XtX = X'*X;
    T* XtX = new T[d*d];
    dot(X.getData(), X.getData(), XtX, n, d, n, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);


//    [Q,L,U] = tygert_svd(XtX,d);
//    Q = double(Q);
//    L = double(diag(L));
    T *Q = new T[d*d];
    T *L = new T[d];
    T *V = NULL;
    random_svd(XtX, d, d, Q, L, V, d);

    delete[] XtX;

//	Xty = X'*y;
    T* Xty = new T[d*Yd];
    dot(X.getData(), Y.getData(), Xty, n, d, Yn, Yd, d, Yd, CblasTrans, CblasNoTrans, CblasColMajor);

//    if isfield(opt,'W0')
    if(opt.hasOpt("W0"))
    {
//        Xty = Xty + opt.W0;
        const gMat2D<T>& W0 = OptMatrix< gMat2D<T> >::dynacast(opt.getOpt("W0"))->getValue();

        if(W0.rows() == d && W0.cols() == Yd)
            axpy(d*Yd, (T)1.0, W0.getData(), 1, Xty, 1);
    }


//    cfr.W = rls_eigen(Q, L, Xty, lambda,d);
    gMat2D<T>* W = new gMat2D<T>(d,Yd);
    T* work = new T[d*(d+1)];
    rls_eigen(Q, L, Xty, W->getData(), lambda, d, d, d, d, d, Yd, work);

    delete [] work;
    delete [] Xty;
    delete [] Q;
    delete [] L;

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

//    cfr.C = [];
    gMat2D<T>* emptyC = new gMat2D<T>();
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

    //	cfr.X = [];
    gMat2D<T>* emptyX = new gMat2D<T>();
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));

    return optimizer;
}


}
#endif // _GURLS_RLSPRIMALR_H_
