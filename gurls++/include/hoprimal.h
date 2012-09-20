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
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);

protected:
    /**
     * Auxiliary method used to call the right eig/svd function for this class
     */
    virtual void eig_function(T* A, T* L, int A_rows_cols, unsigned long d, const GurlsOptionsList &opt);
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
    virtual void eig_function(T* A, T* L, int A_rows_cols, unsigned long d, const GurlsOptionsList &opt);
};





template<typename T>
void ParamSelHoPrimal<T>::eig_function(T* A, T* L, int A_rows_cols,unsigned long , const GurlsOptionsList &)
{
    eig_sm(A, L, A_rows_cols);
}

template<typename T>
void ParamSelHoPrimalr<T>::eig_function(T* A, T* L, int A_rows_cols, unsigned long d, const GurlsOptionsList &opt)
{
    T* V = NULL;
    T k = gurls::round((opt.getOptAsNumber("eig_percentage")*d)/100.0);
    random_svd(A, A_rows_cols, A_rows_cols, A, L, V, A_rows_cols, k);
}

template <typename T>
GurlsOptionsList *ParamSelHoPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //    [n,T]  = size(y);
    const unsigned long n = Y.rows();
    const unsigned long t = Y.cols();

    const unsigned long x_rows = X.rows();
    const unsigned long d = X.cols();


    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");
//    nestedOpt->copyOpt<T>("kernel", opt);
//    nestedOpt->addOpt("predkernel", new GurlsOptionsList("predkernel"));


    const GurlsOptionsList* split = opt.getOptAs<GurlsOptionsList>("split");
    const gMat2D< unsigned long > &indices_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("indices");
    const gMat2D< unsigned long > &lasts_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("lasts");

    const unsigned long *lasts = lasts_mat.getData();
    const unsigned long* indices_buffer = indices_mat.getData();


    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));

    int nholdouts = static_cast<int>(std::ceil( opt.getOptAsNumber("nholdouts")));

    gMat2D<T> *LAMBDA = new gMat2D<T>(1, t);
    T* lambdas = LAMBDA->getData();
    set(lambdas, (T)0.0, t);


    T* Q = new T[d*d];
    T* Qty = new T[d*t];
    T *L = new T[d];

    gMat2D<T>* perf_mat = new gMat2D<T>(nholdouts, tot*t);
    T* perf = perf_mat->getData();

    gMat2D<T>*  guesses_mat = new gMat2D<T>(nholdouts, tot);
    T* ret_guesses = guesses_mat->getData();

    gMat2D<T>* lambdas_round_mat = new gMat2D<T>(nholdouts, t);
    T* lambdas_round = lambdas_round_mat->getData();

    PredPrimal< T > primal;
    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");
    nestedOpt->addOpt("optimizer",optimizer);

    gMat2D<T> *W = new gMat2D<T>(d, t);
    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    for(int nh=0; nh<nholdouts; ++nh)
    {
        unsigned long last = lasts[nh];
        unsigned long* tr = new unsigned long[last];
        unsigned long* va = new unsigned long[n-last];

        //copy int tr indices_ from n*nh to last
        copy< unsigned long >(tr,indices_buffer + n*nh,last,1,1);

        //copy int va indices_ from n*nh+last to n*nh+n
        copy< unsigned long >(va,(indices_buffer+ n*nh+last), n-last,1,1);

        //       K = X(tr,:)'*X(tr,:);
        T* Xtr = new T[last*d];
        subMatrixFromRows(X.getData(), x_rows, d, tr, last, Xtr);

        dot(Xtr, Xtr, Q, last, d, last, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);

        eig_function(Q, L, d, d, opt);

        unsigned long miN = std::min(d,last);
        T* guesses = lambdaguesses(L, d, miN,last, tot, (T)(opt.getOptAsNumber("smallnumber")));

        T* ap = new T[tot*t];

        T* Ytr = new T[last*t];
        subMatrixFromRows(Y.getData(), n, t, tr, last, Ytr);


        T* QtyT = new T[d * last];
        dot(Q, Xtr, QtyT, d, d, last, d, d, last, CblasTrans, CblasTrans, CblasColMajor);

        delete [] Xtr;

        dot(QtyT, Ytr, Qty, d, last, last, t, d, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

        delete [] Ytr;
        delete [] QtyT;


        gMat2D<T> xx(n-last, d);
        gMat2D<T> yy(n-last, t);

        subMatrixFromRows(X.getData(), x_rows, d, va, n-last, xx.getData());
        subMatrixFromRows(Y.getData(), n, t, va, n-last, yy.getData());

        T* work = new T[d*(d+1)];

        for(int i=0; i<tot; ++i)
        {
            rls_eigen(Q, L, Qty, W->getData(), guesses[i], last, d, d, d, d, t, work);

            OptMatrix<gMat2D<T> > *ret_pred = primal.execute(xx, yy, *nestedOpt);

            nestedOpt->removeOpt("pred");
            nestedOpt->addOpt("pred", ret_pred);

            GurlsOptionsList* ret_perf = perfClass->execute(xx, yy, *nestedOpt);

            gMat2D<T> &forho_vec = ret_perf->getOptValue<OptMatrix<gMat2D<T> > >("forho");

            copy(ap+i, forho_vec.getData(), t, tot, 1);

            delete ret_perf;
        }

        delete [] va;
        delete [] tr;
        delete [] work;

        //[dummy,idx] = max(ap,[],1);
        work = NULL;
        unsigned long* idx = new unsigned long[t];
        indicesOfMax(ap, tot, t, idx, work, 1);

        //vout.lambdas_round{nh} = guesses(idx);
        T* lambdas_nh = new T[t];
        copyLocations(idx, guesses, t, tot, lambdas_nh);

        copy(lambdas_round+nh, lambdas_nh, t, nholdouts, 1);

        //add lambdas_nh to lambdas
        axpy(t, (T)1, lambdas_nh, 1, lambdas, 1);

        delete [] lambdas_nh;
        delete [] idx;

        //  vout.perf{nh} = ap;
        copy(perf + nh, ap, tot*t, nholdouts, 1);

        //  vout.guesses{nh} = guesses;
        copy(ret_guesses + nh, guesses, tot, nholdouts, 1);

        delete [] guesses;
        delete [] ap;

    }

    delete nestedOpt;

    delete perfClass;
    delete [] Q;
    delete [] Qty;
    delete [] L;

    GurlsOptionsList* paramsel;

    if(opt.hasOpt("paramsel"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt("paramsel", opt);

        paramsel = GurlsOptionsList::dynacast(tmp_opt->getOpt("paramsel"));
        tmp_opt->removeOpt("paramsel", false);
        delete tmp_opt;

        paramsel->removeOpt("guesses");
        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("perf");
        paramsel->removeOpt("lambdas_round");
    }
    else
        paramsel = new GurlsOptionsList("paramsel");


    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

    if(nholdouts>1)
        scal(t, (T)1.0/nholdouts, lambdas, 1);

    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*LAMBDA));
    paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*perf_mat));
    paramsel->addOpt("lambdas_round", new OptMatrix<gMat2D<T> >(*lambdas_round_mat));

    return paramsel;
}

}

#endif // _GURLS_HOPRIMAL_H_
