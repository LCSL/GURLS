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


#ifndef _GURLS_RLSDUAL_H_
#define _GURLS_RLSDUAL_H_

#include "optimization.h"

#include <set>

namespace gurls {

/**
 * \ingroup Optimization
 * \brief RLSDual is the sub-class of Optimizer that implements RLS with the dual formulation
 */

template <typename T>
class RLSDual: public Optimizer<T>{

public:
    /**
     * Computes a classifier for the dual formulation of RLS.
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
GurlsOptionsList* RLSDual<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)
{
   //	lambda = opt.singlelambda(opt.paramsel.lambdas);
   const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
   T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());

   const GurlsOptionsList* kernel = opt.getOptAs<GurlsOptionsList>("kernel");
   const gMat2D<T>& K_mat = kernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

   T* K = new T[K_mat.getSize()];
   copy(K, K_mat.getData(), K_mat.getSize());

    //n = size(opt.kernel.K,1);
   const long n = K_mat.rows();

   //T = size(y,2);
   const long t = Y.cols();


//    std::cout << "Solving dual RLS... " << std::endl;

    const T coeff = n*static_cast<T>(lambda);
    long i=0;
    for(T* it = K; i<n; ++i, it+=n+1)
        *it += coeff;


   std::set<T*> garbage;

   gMat2D<T>* retC = NULL;

   try // Try solving it with cholesky first.
   {
//        R = chol(K);
        T* R = new T[n*n];
        garbage.insert(R);
        cholesky(K, n, n, R);

//        cfr.C = R\(R'\y);
        retC = new gMat2D<T>(Y.rows(), t);

        copy(retC->getData(), Y.getData(), Y.getSize());
        mldivide_squared(R, retC->getData(), n, n, retC->rows(), retC->cols(), CblasTrans);
        mldivide_squared(R, retC->getData(), n, n, retC->rows(), retC->cols(), CblasNoTrans);

        delete[] R;
        garbage.erase(R);
   }
   catch (gException& /*gex*/)
   {
       for(typename std::set<T*>::iterator it = garbage.begin(); it != garbage.end(); ++it)
           delete[] (*it);

       garbage.clear();

       if(retC != NULL)
           delete retC;


//           [Q,L,V] = svd(K);
//           Q = double(Q);
//           L = double(diag(L));
       T *Q, *L, *V;
       int Q_rows, Q_cols;
       int L_len;
       int V_rows, V_cols;

       svd(K, Q, L, V, n, n, Q_rows, Q_cols, L_len, V_rows, V_cols);


//           cfr.C = rls_eigen(Q,L,y,lambda,n);
       retC = new gMat2D<T>(Q_rows, Y.cols());

       T* work = new T[L_len*(Q_rows+1)];
       rls_eigen(Q, L, Y.getData(), retC->getData(), lambda, n, Q_rows, Q_cols, L_len, Y.rows(), Y.cols(), work);

       delete [] work;
       delete [] Q;
       delete [] L;
       delete [] V;
   }

   delete[] K;

   GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

//       if strcmp(opt.kernel.type, 'linear')
   if(kernel->getOptAsString("type") == "linear")
   {
//           cfr.W = X'*cfr.C;
       gMat2D<T>* W  = new gMat2D<T>(X.cols(), retC->cols());
       dot(X.getData(), retC->getData(), W->getData(), X.rows(), X.cols(), retC->rows(), retC->cols(), W->rows(), W->cols(), CblasTrans, CblasNoTrans, CblasColMajor);
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
#endif // _GURLS_RLSDUAL_H_
