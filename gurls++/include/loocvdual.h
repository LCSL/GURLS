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


#ifndef _GURLS_LOOCVDUAL_H_
#define _GURLS_LOOCVDUAL_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <set>

#include "options.h"
#include "optlist.h"
#include "gmat2d.h"
#include "gvec.h"
#include "gmath.h"

#include "paramsel.h"
#include "perf.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelLoocvDual is the sub-class of ParamSelection that implements LOO cross-validation with the dual formulation
 */

template <typename T>
class ParamSelLoocvDual: public ParamSelection<T>{

public:
    /**
     * Performs parameter selection when the dual formulation of RLS is used.
     * The leave-one-out approach is used.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - nlambda (default)
     *  - hoperf (default)
     *  - smallnumber (default)
     *  - kernel (settable with the class Kernel and its subclasses)
     *
     * \return adds the following fields to opt:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - guesses = array of guesses for the regularization parameter lambda
     *  - acc = matrix of validation accuracies for each lambda guess and for each class
     */
   void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void ParamSelLoocvDual<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{
//    [n,T]  = size(y);
    const int n = Y_OMR.rows();
    const int t = Y_OMR.cols();
    const int d = X_OMR.cols();

    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);


//    tot = opt.nlambda;
    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));


//    [Q,L] = eig(opt.kernel.K);
//    Q = double(Q);
//    L = double(diag(L));
    GurlsOptionsList* kernel = GurlsOptionsList::dynacast(opt.getOpt("kernel"));
    GurlsOption *K_opt = kernel->getOpt("K");


    gMat2D<T> *K_mat = &(OptMatrix<gMat2D<T> >::dynacast(K_opt))->getValue();
    gMat2D<T> K(K_mat->cols(), K_mat->rows());
    K_mat->transpose(K);

    const int qrows = K_mat->rows();
    const int qcols = K_mat->cols();
    const int l_length = qrows;

    T *Q = K.getData();
    T *L = new T[l_length];

    eig_sm(Q, L, qrows); // qrows == qcols

    int r = n;
    if(kernel->getOptAsString("type") == "linear")
    {
      //set(L+X_OMR.cols(), (T) 1.0e-12, l_length-X_OMR.cols());
        set(L, (T) 1.0e-12, l_length-X_OMR.cols());
        r = std::min(n,d);
    }


//    Qty = Q'*y;
//    T* Qty = new T[qrows*qcols];
    T* Qty = new T[qcols*t];
    dot(Q, Y.getData(), Qty, qrows, qcols, n, t, qcols, t, CblasTrans, CblasNoTrans, CblasColMajor);

    T* guesses = lambdaguesses(L, n, r, n, tot, (T)(opt.getOptAsNumber("smallnumber")));


    GurlsOption* pred_old = NULL;
    GurlsOption* perf_old = NULL;

    if(opt.hasOpt("pred"))
    {
        pred_old = opt.getOpt("pred");
        opt.removeOpt("pred", false);
    }

    gMat2D<T>* pred = new gMat2D<T>(n, t);
    OptMatrix<gMat2D<T> >* pred_opt = new OptMatrix<gMat2D<T> >(*pred);
    opt.addOpt("pred", pred_opt);

    if(opt.hasOpt("perf"))
    {
        perf_old = opt.getOpt("perf");
        opt.removeOpt("perf", false);
    }

    GurlsOptionsList* perf = new GurlsOptionsList("perf");
    opt.addOpt("perf", perf);


    gMat2D<T> tmp_pred(pred->cols(), pred->rows());


    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));

    T* ap = new T[tot*t];
    T* C_div_Z = new T[qrows];
    T* C = new T[qrows*qcols];
    T* Z = new T[qrows];

    for(int i = 0; i < tot; ++i)
    {
        rls_eigen(Q, L, Qty, C, guesses[i], n, qrows, qcols, l_length, qcols, t);
        GInverseDiagonal(Q, L, guesses+i, Z, qrows, qcols, l_length, 1);

        for(int j = 0; j< t; ++j)
        {
            rdivide(C + (qrows*j), Z, C_div_Z, qrows);

//            opt.pred(:,t) = y(:,t) - (C(:,t)./Z);
            copy(tmp_pred.getData()+(n*j), Y.getData() + (n*j), n);
            axpy(n, (T)-1.0, C_div_Z, 1, tmp_pred.getData() + (n*j), 1);
        }

        tmp_pred.transpose(*pred);

//        opt.perf = opt.hoperf([],y,opt);
        const gMat2D<T> dummy;
        perfClass->execute(dummy, Y_OMR, opt);

        gMat2D<T> *forho_vec = &(OptMatrix<gMat2D<T> >::dynacast(perf->getOpt("forho")))->getValue();

        copy(ap+i, forho_vec->getData(), t, tot, 1);

    }

    delete [] C;
    delete [] Z;
    delete [] C_div_Z;
    delete perfClass;

    delete[] L;
    //delete[] Q;

    unsigned long* idx = new unsigned long[t];
    T* work = NULL;
    indicesOfMax(ap, tot, t, idx, work, 1);

    T* lambdas = copyLocations(idx, guesses, t, tot);

    delete[] idx;

    OptNumberList* LAMBDA = new OptNumberList();
    for (T* l_it = lambdas, *l_end = lambdas+t; l_it != l_end; ++l_it)
        LAMBDA->add(static_cast<double>(*l_it));

    delete[] lambdas;

//    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");
    GurlsOptionsList* paramsel = NULL;
    if(!opt.hasOpt("paramsel"))
    {
        paramsel = new GurlsOptionsList("paramsel");
        opt.addOpt("paramsel", paramsel);
    }
    else
    {
        paramsel = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));

        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("perf");
        paramsel->removeOpt("guesses");

    }

//     opt.addOpt("lambdas", LAMBDA);
    paramsel->addOpt("lambdas", LAMBDA);

    //vout.perf = 	ap;

    gMat2D<T>* looe_mat = new gMat2D<T>(tot, t);
    transpose(ap, tot, t, looe_mat->getData());

    delete[] ap;

    paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*looe_mat));

    //vout.guesses = 	guesses;
    gMat2D<T> *guesses_mat = new gMat2D<T>(guesses, 1, tot, true);
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

    delete[] guesses;

    opt.removeOpt("pred");
    if(pred_old != NULL)
        opt.addOpt("pred", pred_old);

    opt.removeOpt("perf");
    if(perf_old != NULL)
        opt.addOpt("perf", perf_old);
}


}

#endif // _GURLS_LOOCVDUAL_H_
