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


#ifndef _GURLS_HOPRIMAL_H_
#define _GURLS_HOPRIMAL_H_

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
#include "dual.h"
#include "utils.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelHoPrimal is the subclass of ParamSelection that implements hold-out cross validation with the primal formulation of RLS
 */

template <typename T>
class ParamSelHoPrimal: public ParamSelection<T>{

public:
    /**
     * Performs parameter selection when the primal formulation of RLS is used.
     * The hold-out approach is used.
     * The performance measure specified by opt.hoperf is maximized.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - nlambda (default)
     *  - hoperf (default)
     *  - smallnumber (default)
     *  - split (settable with the class Split and its subclasses)
     *
     * \return adds the following fields to opt:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - guesses = array of guesses for the regularization parameter lambda
     *  - forho = matrix of validation accuracies for each lambda guess and for each class
     */
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);

protected:
    /**
     * Auxiliary method used to call the right eig/svd function for this class
     */
    virtual void eig_function(T* A, T* L, int A_rows_cols);
};


/**
 * \ingroup ParameterSelection
 * \brief ParamSelHoPrimalr is the randomized version of \ref ParamSelHoPrimal
 */

template <typename T>
class ParamSelHoPrimalr: public ParamSelHoPrimal<T>{

protected:
    /**
     * Auxiliary method used to call the right eig/svd function for this class
     */
    virtual void eig_function(T* A, T* L, int A_rows_cols);
};





template<typename T>
void ParamSelHoPrimal<T>::eig_function(T* A, T* L, int A_rows_cols)
{
    eig_sm(A, L, A_rows_cols);
}

template<typename T>
void ParamSelHoPrimalr<T>::eig_function(T* A, T* L, int A_rows_cols)
{
    T* V = NULL;
    random_svd(A, A_rows_cols, A_rows_cols, A, L, V, A_rows_cols);
}

template <typename T>
void ParamSelHoPrimal<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{
    //    [n,T]  = size(y);
    const unsigned long n = Y_OMR.rows();
    const unsigned long t = Y_OMR.cols();

    const unsigned long x_rows = X_OMR.rows();
    const unsigned long d = X_OMR.cols();

    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);

    GurlsOptionsList* perf_old = NULL;
    if(opt.hasOpt("perf"))
    {
        perf_old = GurlsOptionsList::dynacast(opt.getOpt("perf"));
        opt.removeOpt("perf", false);
    }
    GurlsOption* pred_old = NULL;

    if(opt.hasOpt("pred"))
    {
        pred_old = opt.getOpt("pred");
        opt.removeOpt("pred", false);
    }
    GurlsOptionsList* optimizer_old = NULL;
    if(opt.hasOpt("optimizer"))
    {
        optimizer_old = GurlsOptionsList::dynacast(opt.getOpt("optimizer"));
        opt.removeOpt("optimizer", false);
    }

    GurlsOptionsList* split = GurlsOptionsList::dynacast(opt.getOpt("split"));
    GurlsOption *index_opt = split->getOpt("indices");
    GurlsOption *lasts_opt = split->getOpt("lasts");


    gMat2D<unsigned long > *indices_mat = &(OptMatrix<gMat2D< unsigned long > >::dynacast(index_opt))->getValue();
    //vector of int one for column of indices_mat
    gMat2D< unsigned long > *lasts_mat = &(OptMatrix<gMat2D< unsigned long > >::dynacast(lasts_opt))->getValue();
    unsigned long *lasts_ = lasts_mat->getData();

    gMat2D< unsigned long > Indices(indices_mat->cols(), indices_mat->rows());
    indices_mat->transpose(Indices);
    unsigned long* indices_buffer = new unsigned long[indices_mat->cols()*indices_mat->rows()];
    copy< unsigned long >(indices_buffer,Indices.getData(),indices_mat->cols()*indices_mat->rows());

    unsigned long *indices_from0toD = new unsigned long[d];
    unsigned long *indices_from0toT = new unsigned long[t];

    for(unsigned long i=0; i<d; ++i)
        indices_from0toD[i]=i;
    for(unsigned long i=0;i<t;++i)
        indices_from0toT[i]=i;

    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));

    int nholdouts = static_cast<int>(std::ceil( opt.getOptAsNumber("nholdouts")));
    T* lambdas = new T[t];
    set(lambdas, (T)0.0, t);


    T* Q = new T[d*d];
    T* Qty = new T[d*t];
    T *L = new T[d];
    T* W = new T[d*t];

    gMat2D<T>* perf_mat = new gMat2D<T>(nholdouts, tot*t);
    T* perf = perf_mat->getData();

    gMat2D<T>*  guesses_mat = new gMat2D<T>(nholdouts, tot);
    T* ret_guesses = guesses_mat->getData();

    gMat2D<T>* lambdas_round_mat = new gMat2D<T>(nholdouts, t);
    T* lambdas_round = lambdas_round_mat->getData();

    PredPrimal< T > primal;
    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));


    for(int nh=0; nh<nholdouts; ++nh)
    {
        unsigned long lasts = lasts_[nh];
        unsigned long* tr = new unsigned long[lasts];
        unsigned long* va = new unsigned long[n-lasts];

        //copy int tr indices_ from n*nh to lasts
        copy< unsigned long >(tr,indices_buffer + n*nh,lasts,1,1);

        //copy int va indices_ from n*nh+lasts to n*nh+n
        copy< unsigned long >(va,(indices_buffer+ n*nh+lasts), n-lasts,1,1);

        //       K = X(tr,:)'*X(tr,:);
        T* Xtr = new T[lasts*d];
        copy_submatrix(Xtr, X.getData(), x_rows, lasts, d, tr, indices_from0toD);

        //Get K(tr,tr) from k_buffer_total
        dot(Xtr, Xtr, Q, lasts, d, lasts, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);

        eig_function(Q, L, d);

        unsigned long miN = std::min(d,lasts);
        T* guesses = lambdaguesses(L, d, miN,lasts, tot, (T)(opt.getOptAsNumber("smallnumber")));

        T* ap = new T[tot*t];

        T* Ytr = new T[lasts*t];
        copy_submatrix(Ytr, Y.getData(), n, lasts, t, tr, indices_from0toT);


        T* QtyT = new T[d * lasts];
        dot(Q, Xtr, QtyT, d, d, lasts, d, d, lasts, CblasTrans, CblasTrans, CblasColMajor);

        dot(QtyT, Ytr, Qty, d, lasts, lasts, t, d, t, CblasNoTrans, CblasNoTrans, CblasColMajor);
        delete [] QtyT;

        GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");
        opt.addOpt("optimizer",optimizer);

        for(int i=0; i<tot; ++i)
        {
            rls_eigen(Q, L, Qty, W, guesses[i], lasts, d, d, d, d, t);

            gMat2D<T> w_tmp(W, t, d, false);
            gMat2D<T>* w_t = new gMat2D<T>(d, t);
            w_tmp.transpose(*w_t);

            optimizer->removeOpt("W");
            optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*w_t));

            T* Xva = new T[(n-lasts)*d];
            copy_submatrix(Xva, X.getData(), x_rows, n-lasts, d, va, indices_from0toD);

            T* Yva = new T[(n-lasts)*t];
            copy_submatrix(Yva, Y.getData(), n, n-lasts, t, va, indices_from0toT);

            const gMat2D<T> Xva_m(Xva,d, n-lasts, false);
            gMat2D<T> xx(n-lasts, d);
            Xva_m.transpose(xx);

            const gMat2D<T> Yva_m(Yva, t, n-lasts, false);
            gMat2D<T> yy(n-lasts, t);
            Yva_m.transpose(yy);

            primal.execute(xx, yy, opt);

            perfClass->execute(xx, yy, opt);

            delete [] Yva;
            delete [] Xva;

            GurlsOptionsList* perf_opt = GurlsOptionsList::dynacast(opt.getOpt("perf"));
            gMat2D<T> &forho_vec = OptMatrix<gMat2D<T> >::dynacast(perf_opt->getOpt("forho"))->getValue();

            copy(ap+i, forho_vec.getData(), t, tot, 1);
        }

        delete perfClass;
        delete [] va;
        delete [] tr;
        delete [] Xtr;
        delete [] Ytr;

        //[dummy,idx] = max(ap,[],1);
        T* work = NULL;
        unsigned long* idx = new unsigned long[t];
        indicesOfMax(ap, tot, t, idx, work, 1);

        //vout.lambdas_round{nh} = guesses(idx);
        T* lambdas_nh = copyLocations(idx, guesses, t, tot);

        copy(lambdas_round+nh*t, lambdas_nh, t);

        //add lambdas_nh to lambdas
        axpy(t, (T)1, lambdas_nh, 1, lambdas, 1);

        delete [] lambdas_nh;
        delete [] idx;

        //  vout.perf{nh} = ap;
        copy(perf + nh*tot*t, ap, tot*t);

        //  vout.guesses{nh} = guesses;
        copy(ret_guesses + nh*tot, guesses, tot);

        delete [] guesses;
        delete [] ap;

    }

    delete [] indices_buffer;
    delete [] Q;
    delete [] Qty;
    delete [] L;
    delete [] W;


    GurlsOptionsList* paramsel = NULL;
    if(!opt.hasOpt("paramsel"))
    {
        paramsel = new GurlsOptionsList("paramsel");
        opt.addOpt("paramsel", paramsel);
    }
    else
    {
        paramsel = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));

        paramsel->removeOpt("guesses");
        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("perf");
        paramsel->removeOpt("lambdas_round");
    }


    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

    if(nholdouts>1)
        scal(t, (T)1.0/nholdouts, lambdas, 1);

    OptNumberList* LAMBDA = new OptNumberList();
    for (T* l_it = lambdas, *l_end = lambdas+t; l_it != l_end; ++l_it)
        LAMBDA->add(static_cast<double>(*l_it));

    delete [] lambdas;

    paramsel->addOpt("lambdas", LAMBDA);
    paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*perf_mat));
    paramsel->addOpt("lambdas_round", new OptMatrix<gMat2D<T> >(*lambdas_round_mat));

    opt.removeOpt("pred");
    if(pred_old != NULL)
        opt.addOpt("pred", pred_old);

    opt.removeOpt("perf");
    if(perf_old != NULL)
        opt.addOpt("perf", perf_old);

    opt.removeOpt("optimizer");
    if(optimizer_old != NULL)
        opt.addOpt("optimizer", optimizer_old);

}

}

#endif // _GURLS_HOPRIMAL_H_
