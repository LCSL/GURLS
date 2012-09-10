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
   GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList* ParamSelLoocvDual<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
//    [n,T]  = size(y);
    const unsigned long n = Y.rows();
    const unsigned long t = Y.cols();

    const unsigned long d = X.cols();


//    tot = opt.nlambda;
    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));


//    [Q,L] = eig(opt.kernel.K);
//    Q = double(Q);
//    L = double(diag(L));
    const GurlsOptionsList* kernel = opt.getOptAs<GurlsOptionsList>("kernel");

    const gMat2D<T> &K_mat = kernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

    gMat2D<T> K(K_mat.rows(), K_mat.cols());
    copy(K.getData(), K_mat.getData(), K_mat.getSize());

    const unsigned long qrows = K.rows();
    const unsigned long qcols = K.cols();
    const unsigned long l_length = qrows;

    T *Q = K.getData();
    T *L = new T[l_length];

    eig_sm(Q, L, qrows); // qrows == qcols

    int r = n;
    if(kernel->getOptAsString("type") == "linear")
    {
        set(L, (T) 1.0e-12, l_length-d);
        r = std::min(n,d);
    }


//    Qty = Q'*y;
//    T* Qty = new T[qrows*qcols];
    T* Qty = new T[qcols*t];
    dot(Q, Y.getData(), Qty, qrows, qcols, n, t, qcols, t, CblasTrans, CblasNoTrans, CblasColMajor);

    T* guesses = lambdaguesses(L, n, r, n, tot, (T)(opt.getOptAsNumber("smallnumber")));



    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");

    gMat2D<T>* pred = new gMat2D<T>(n, t);
    OptMatrix<gMat2D<T> >* pred_opt = new OptMatrix<gMat2D<T> >(*pred);
    nestedOpt->addOpt("pred", pred_opt);


    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));

    gMat2D<T>* perf = new gMat2D<T>(tot, t);
    T* ap = perf->getData();

    T* C_div_Z = new T[qrows];
    T* C = new T[qrows*qcols];
    T* Z = new T[qrows];

    for(int i = 0; i < tot; ++i)
    {
        rls_eigen(Q, L, Qty, C, guesses[i], n, qrows, qcols, l_length, qcols, t);
        GInverseDiagonal(Q, L, guesses+i, Z, qrows, qcols, l_length, 1);

        for(unsigned long j = 0; j< t; ++j)
        {
            rdivide(C + (qrows*j), Z, C_div_Z, qrows);

//            opt.pred(:,t) = y(:,t) - (C(:,t)./Z);
            copy(pred->getData()+(n*j), Y.getData() + (n*j), n);
            axpy(n, (T)-1.0, C_div_Z, 1, pred->getData() + (n*j), 1);
        }

//        opt.perf = opt.hoperf([],y,opt);
        const gMat2D<T> dummy;
        GurlsOptionsList* perf = perfClass->execute(dummy, Y, *nestedOpt);

        gMat2D<T> &forho_vec = perf->getOptValue<OptMatrix<gMat2D<T> > >("forho");

        copy(ap+i, forho_vec.getData(), t, tot, 1);

        delete perf;
    }

    delete nestedOpt;
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


    gMat2D<T> *LAMBDA = new gMat2D<T>(1, t);
    copy(LAMBDA->getData(), lambdas, t);

    delete[] lambdas;


    GurlsOptionsList* paramsel;

    if(opt.hasOpt("paramsel"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt<T>("paramsel", opt);

        paramsel = GurlsOptionsList::dynacast(tmp_opt->getOpt("paramsel"));
        tmp_opt->removeOpt("paramsel", false);
        delete tmp_opt;

        paramsel->removeOpt("guesses");
        paramsel->removeOpt("perf");
        paramsel->removeOpt("lambdas");
    }
    else
        paramsel = new GurlsOptionsList("paramsel");


//     opt.addOpt("lambdas", LAMBDA);
    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*LAMBDA));

    //vout.perf = 	ap;
    paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*perf));

    //vout.guesses = 	guesses;
    gMat2D<T> *guesses_mat = new gMat2D<T>(guesses, 1, tot, true);
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

    delete[] guesses;

    return paramsel;

}


}

#endif // _GURLS_LOOCVDUAL_H_
