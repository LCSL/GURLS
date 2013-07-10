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


#ifndef _GURLS_RLSDUALR_H_
#define _GURLS_RLSDUALR_H_

#include "gurls++/optimization.h"
#include "gurls++/utils.h"
#include "gurls++/gmath.h"

#include <set>

namespace gurls {

/**
 * \ingroup Optimization
 * \brief RLSDualr is the sub-class of Optimizer that implements RLS with the dual formulation, using a randomized version of SVD
 */

template <typename T>
class RLSDualr: public Optimizer<T>{

public:
    /**
     * Computes a classifier for the dual formulation of RLS, using a randomized version of Singular value decomposition.
     * The regularization parameter is set to the one found in the field paramsel of opt.
     * In case of multiclass problems, the regularizers need to be combined with the function specified inthe field singlelambda of opt
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - singlelambda (default)
     *  - paramsel (settable with the class ParamSelection and its subclasses)
     *  - kernel (settable with the class Kernel and its subclasses)
     *
     * \return adds to opt the field optimizer, which is a list containing the following fields:
     *  - W = empty matrix
     *  - C = matrix of coefficient vectors of dual rls estimator for each class
     *  - X = empty matrix
     */
   GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* RLSDualr<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)
{
//	lambda = opt.singlelambda(opt.paramsel.lambdas);
    const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
    T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());

    const GurlsOptionsList* kernel = opt.getOptAs<GurlsOptionsList>("kernel");
    const gMat2D<T>& K_mat = kernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

    T* K = new T[K_mat.getSize()];
    copy(K, K_mat.getData(), K_mat.getSize());

    //n = size(opt.kernel.K,1);
    const unsigned long n = K_mat.rows();

   //T = size(y,2);
    const unsigned long t = Y.cols();

    const T coeff = n*lambda;
    unsigned long i=0;
    for(T* it = K; i<n; ++i, it+=n+1)
        *it += coeff;


//       [Q,L,V] = tygert_svd(K,n);
//       Q = double(Q);
//       L = double(diag(L));
    T *Q = new T[n*n];
    T *L = new T[n];
    T *V = NULL;

//    k = round(opt.eig_percentage*n/100);
    unsigned long k = static_cast<unsigned long>(gurls::round((opt.getOptAsNumber("eig_percentage")*n)/100.0));
    random_svd(K, n, n, Q, L, V, k);


    gMat2D<T> *retC = new gMat2D<T>(n,t);
    T* work = new T[n*(n+1)];

    T* Qty = new T[n*Y.cols()];
    dot(Q, Y.getData(), Qty, n, n, Y.rows(), Y.cols(), n, Y.cols(), CblasTrans, CblasNoTrans, CblasColMajor);

    rls_eigen(Q, L, Qty, retC->getData(), lambda, n, n, n, n, n, t, work);

    delete [] Qty;
    delete [] work;
    delete [] Q;
    delete [] L;

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");


//       if strcmp(opt.kernel.type, 'linear')
   if(kernel->getOptAsString("type") == "linear")
   {
   //           cfr.W = X'*cfr.C;
        gMat2D<T>* W  = new gMat2D<T>(X.cols(), t);
        dot(X.getData(), retC->getData(), W->getData(), X.rows(), X.cols(), n, t, X.cols(), t, CblasTrans, CblasNoTrans, CblasColMajor);

        optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

//           cfr.C = [];
        gMat2D<T>* emptyC = new gMat2D<T>();
        optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

//           cfr.X = [];
        gMat2D<T>* emptyX = new gMat2D<T>();
        optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));

        delete retC;
    }
    else
    {
//           cfr.W = [];
        gMat2D<T>* emptyW = new gMat2D<T>();
        optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*emptyW));

//           cfr.C = retC;
        optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*retC));

//           cfr.X = X;
        gMat2D<T>* optX = new gMat2D<T>(X);
        optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*optX));
    }

    return optimizer;
}

}
#endif // _GURLS_RLSDUALR_H_
