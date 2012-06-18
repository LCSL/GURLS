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
   void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};


template <typename T>
void RLSDual<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{
   //	lambda = opt.singlelambda(opt.paramsel.lambdas);
//    std::vector<double> ll = OptNumberList::dynacast(opt.getOpt("lambdas"))->getValue();
   GurlsOptionsList* paramsel = static_cast<GurlsOptionsList*>(opt.getOpt("paramsel"));
   std::vector<double> ll = OptNumberList::dynacast(paramsel->getOpt("lambdas"))->getValue();
   T lambda = static_cast<T>((OptFunction::dynacast(opt.getOpt("singlelambda")))->getValue(&(*(ll.begin())), ll.size()));

   gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
   X_OMR.transpose(X);

   gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
   Y_OMR.transpose(Y);


   GurlsOptionsList* kernel = static_cast<GurlsOptionsList*>(opt.getOpt("kernel"));
   GurlsOption *K_opt = kernel->getOpt("K");

//    if (K_opt->getDataID() != typeid(T))
//        return;

   gMat2D<T> *K_mat = &(OptMatrix<gMat2D<T> >::dynacast(K_opt))->getValue();
   gMat2D<T> K(K_mat->cols(), K_mat->rows());
   K_mat->transpose(K);

    //n = size(opt.kernel.K,1);
   const long n = K_mat->rows();

   //T = size(y,2);
   const long t = Y_OMR.cols();


    std::cout << "Solving dual RLS... " << std::endl;

    const T coeff = n*static_cast<T>(lambda);
    long i=0;
    for(T* it = K.getData(); i<n; ++i, it+=n+1)
        *it += coeff;


   std::set<T*> garbage;

//    bool retX = false;
   T* retC = NULL;
//   T* retX = NULL;
   int retC_rows = 0;
   int retC_cols = 0;

   try{ // Try solving it with cholesky first.

//		R = chol(K);
       T* R = cholesky(K.getData(), n, n);
       garbage.insert(R);


////           cfr.C = R\(R'\y);
//       T* Rt = transpose_cm(R, n, n);
//       garbage.insert(Rt);


//       int RtMlDivY_rows, RtMlDivY_cols;
//       T* RtMlDivY = mldivide(Rt, Y.getData(), n, n, Y_OMR.rows(), t, RtMlDivY_rows, RtMlDivY_cols);
//       garbage.insert(RtMlDivY);

//       delete[] Rt;
//       garbage.erase(Rt);

//       retC = mldivide(R, RtMlDivY, n, n, RtMlDivY_rows, RtMlDivY_cols, retC_rows, retC_cols);
//       garbage.insert(retC);


//       delete[] R;
//       garbage.erase(R);

//       delete[] RtMlDivY;
//       garbage.erase(RtMlDivY);


       retC_rows = Y_OMR.rows();
       retC_cols = t;

        retC = new T[retC_rows*retC_cols];
        garbage.insert(retC);

        copy(retC, Y.getData(), Y.getSize());
        mldivide_squared(R, retC, n, n, retC_rows, retC_cols, CblasTrans);
        mldivide_squared(R, retC, n, n, retC_rows, retC_cols, CblasNoTrans);


        delete[] R;
        garbage.erase(R);

//           cfr.X = X;
//        retX = true;
//       retX = new T[X_OMR.rows()*X_OMR.cols()];
//       garbage.insert(retX);
//       copy(retX, X.getData(), X.getSize());

   }
   catch (gException& /*gex*/)
   {
       for(typename std::set<T*>::iterator it = garbage.begin(); it != garbage.end(); ++it)
           delete[] (*it);

       garbage.clear();


//           [Q,L,V] = svd(K);
//           Q = double(Q);
//           L = double(diag(L));
       T *Q, *L, *V;
       int Q_rows, Q_cols;
       int L_len;
       int V_rows, V_cols;

       svd(K.getData(), Q, L, V, n, n, Q_rows, Q_cols, L_len, V_rows, V_cols);


//        const int Q_rows = K_mat->rows();
//        const int Q_cols = K_mat->cols();
//        int L_len = Q_rows;

//        T *Q = K.getData();
//        T *L = new T[L_len];

//        eig_sm(Q, L, Q_rows, Q_cols);


//           cfr.C = rls_eigen(Q,L,y,lambda,n);
//       retC = rls_eigen(Q, L, Y.getData(), &lambda, n, Q_rows, Q_cols, L_len, 1, Y_OMR.rows(), Y_OMR.cols());
       retC_rows = Q_rows;
       retC_cols = Y_OMR.cols();

       retC = new T[retC_rows*retC_cols];
       rls_eigen(Q, L, Y.getData(), retC, lambda, n, Q_rows, Q_cols, L_len, Y_OMR.rows(), Y_OMR.cols());

//       delete [] Q;
       delete [] L;
//       delete [] V;
   }

   GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

//       if strcmp(opt.kernel.type, 'linear')
   if(kernel->getOptAsString("type") == "linear")
   {
   //           cfr.W = X'*cfr.C;
        T* retW  = new T[X_OMR.cols()*retC_cols];
        dot(X.getData(), retC, retW, X_OMR.rows(), X_OMR.cols(), retC_rows, retC_cols, X_OMR.cols(), retC_cols, CblasTrans, CblasNoTrans, CblasColMajor);

        gMat2D<T> tmp(retW, retC_cols, X_OMR.cols(), false);
        gMat2D<T>* W = new gMat2D<T>(X_OMR.cols(), retC_cols);
        tmp.transpose(*W);

        delete[] retW;

    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

//           cfr.C = [];
        gMat2D<T>* emptyC = new gMat2D<T>();
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

//           cfr.X = [];
        gMat2D<T>* emptyX = new gMat2D<T>();
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));

   }
   else
   {
//           cfr.W = [];
       gMat2D<T>* emptyW = new gMat2D<T>();
        optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*emptyW));

//           cfr.C = retC;
        gMat2D<T> tmp(retC, retC_cols, retC_rows, false);
        gMat2D<T>* C = new gMat2D<T>(retC_rows, retC_cols);
        tmp.transpose(*C);

        optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*C));

//           cfr.X = X;
//         if(retX)
//         {
            gMat2D<T>* optX = new gMat2D<T>(X_OMR);
            optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*optX));
//         }

   }

    opt.removeOpt("optimizer", false);
    opt.addOpt("optimizer",optimizer);

   delete[] retC;
}

}
#endif // _GURLS_RLSDUAL_H_
