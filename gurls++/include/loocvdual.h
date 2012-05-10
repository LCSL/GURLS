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

template <typename T>
class LoocvDual: public ParamSelection<T>{

public:
   void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void LoocvDual<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{
//    [n,T]  = size(y);
    const int n = Y_OMR.rows();
    const int t = Y_OMR.cols();

    gMat2D<float> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    gMat2D<float> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);


//    tot = opt.nlambda;
    int tot = std::ceil( static_cast<OptNumber*>(opt.getOpt("nlambda"))->getValue() );


//    [Q,L] = eig(opt.kernel.K);
//    Q = double(Q);
//    L = double(diag(L));
    GurlsOptionsList* kernel = static_cast<GurlsOptionsList*>(opt.getOpt("kernel"));
    GurlsOption *K_opt = kernel->getOpt("K");

//    if (K_opt->getDataID() != typeid(T))
//        return;

    gMat2D<T> *K_mat = &(OptMatrix<gMat2D<T> >::dynacast(K_opt))->getValue();
    gMat2D<T> K(K_mat->cols(), K_mat->rows());
    K_mat->transpose(K);

    const int qrows = K_mat->rows();
    const int qcols = K_mat->cols();
    const int l_length = qrows;

    T *Q = K.getData();
    T *L = new T[l_length];

    eig_sm(Q, L, qrows, qcols);

    if(kernel->getOptAsString("type") == "linear")
        //set(L+X_OMR.cols(), (T) 1.0e-12, l_length-X_OMR.cols());
        set(L, (T) 1.0e-12, l_length-X_OMR.cols());


//    Qty = Q'*y;
//    T* Qty = new T[qrows*qcols];
    T* Qty = new T[qcols*t];
    dot(Q, Y.getData(), Qty, qrows, qcols, n, t, qcols, t, CblasTrans, CblasNoTrans, CblasColMajor);

    // WARNING ============================== -> SEE LAMBDAGUESSES
    //	filtered = L(L > 200*eps^0.5);
    // WARNING: using an approximate version of the eps Matlab variable
    //gVec<T> filtered = (L.asMatrix()).where(L.asMatrix() > 2.9802e-6f);
//    T* filtered = L;

    //	lmin = min(filtered)/n;
    //	lmax = max(filtered)/n;
//    const T lmin = (*std::min_element(L, L+qrows))/ n;
//    const T lmax = (*std::max_element(L, L+qrows))/ n;

    //	q = (lmax/lmin)^(1/tot);
//    const T q = pow((lmax/lmin), ( static_cast<T>(1.0) / tot));


//    guesses = zeros(1,tot);
//    T* guesses = new T[tot];
//    set(guesses, 0.f, tot);

    //T* guesses = lambdaguesses(L, n, n, n, tot, (T)(opt.getOptAsNumber("smallnumber")));
    T* guesses = lambdaguesses(L, n, n, n, tot, (T)1.0e8);

//    LOOSQE = zeros(tot,T);
    T* LOOSQE = new T[tot*t];
    set(LOOSQE, 0.f, tot*t);


    GurlsOption* pred_old = NULL;

    if(opt.hasOpt("pred"))
    {
        pred_old = opt.getOpt("pred");
        opt.removeOpt("pred", false);
    }

    gMat2D<T>* pred = new gMat2D<T>(n, t);
    OptMatrix<gMat2D<T> >* pred_opt = new OptMatrix<gMat2D<T> >(*pred);
    opt.addOpt("pred", pred_opt);

    GurlsOptionsList* perf = new GurlsOptionsList("perf");
    opt.addOpt("perf", perf);

//    if (pred_opt->getDataID() != typeid(T))
//        return;

    gMat2D<T> tmp_pred(pred->cols(), pred->rows());

    //const int pred_size = pred->getSize();
    //T* tmp_pred = new T[pred_size];

    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));

    T* ap = new T[tot*t];
    T* C_div_Z = new T[qrows];
    T* C = new T[qrows*qcols];
    T* Z = new T[qrows];

//    for i = 1:tot
    for(int i = 0; i < tot; ++i)
    {
//        guesses(i) = lmin*(q^i);
//        guesses[i] = lmin*std::pow(q,i);

//        C = rls_eigen(Q,L,Qty,guesses(i),n);
        //C.size() == qrowsXqcols
        rls_eigen(Q, L, Qty, C, guesses[i], n, qrows, qcols, l_length, qcols, t);

//        Z = GInverseDiagonal(Q,L,guesses(i));
        //Z.size() == qrowsX1
        GInverseDiagonal(Q, L, guesses+i, Z, qrows, qcols, l_length, 1);

//        opt.pred = zeros(n,T);
        //set(tmp_pred.getData(), (T)0.0, pred_size);

//        for t = 1:T
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

    //[dummy,idx] = max(ap,[],1);
    const unsigned int* idx = indicesOfMax(ap, tot, t, 1);

    //vout.lambdas = 	guesses(idx);
    T* lambdas = copyLocations(idx, guesses, t, tot);

    delete[] idx;

    OptNumberList* LAMBDA = new OptNumberList();
    for (T* l_it = lambdas, *l_end = lambdas+t; l_it != l_end; ++l_it)
        LAMBDA->add(static_cast<double>(*l_it));

    delete[] lambdas;

    opt.addOpt("lambdas", LAMBDA);

    //vout.looe{1} = 	ap;

    T* tmp = transpose_rm(ap, t, tot);

    delete[] ap;

    gMat2D<T>* looe_mat = new gMat2D<T>(tmp, tot, t, true);

    delete[] tmp;

    opt.addOpt("looe", new OptMatrix<gMat2D<T> >(*looe_mat));


    //vout.guesses = 	guesses;
    gVec<T> guessesVector(guesses, tot, false);
    opt.addOpt("guesses", new OptMatrix<gMat2D<T> >(guessesVector.asMatrix()));

    delete[] guesses;

    opt.removeOpt("pred");
    if(pred_old != NULL)
        opt.addOpt("pred", pred_old);
}


}

#endif // _GURLS_LOOCVDUAL_H_
