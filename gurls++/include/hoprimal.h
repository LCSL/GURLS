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
     * \return paramsel, a GurlsOptionList with the following fields:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - guesses = array of guesses for the regularization parameter lambda
     *  - forho = matrix of validation accuracies for each lambda guess and for each class
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);

protected:
    /**
     * Auxiliary method used to call the right eig/svd function for this class
     */
    unsigned long eig_function(T* A, T* L, int A_rows_cols, unsigned long d, const GurlsOptionsList &opt, unsigned long last);
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
    unsigned long eig_function(T* A, T* L, int A_rows_cols, unsigned long d, const GurlsOptionsList &opt, unsigned long last);
};





template<typename T>
unsigned long ParamSelHoPrimal<T>::eig_function(T* A, T* L, int A_rows_cols,unsigned long d, const GurlsOptionsList &, unsigned long last)
{
    eig_sm(A, L, A_rows_cols);

    return std::min(d,last);
}

template<typename T>
unsigned long ParamSelHoPrimalr<T>::eig_function(T* A, T* L, int A_rows_cols, unsigned long d, const GurlsOptionsList &opt, unsigned long )
{
    T* V = NULL;
    T k = gurls::round((opt.getOptAsNumber("eig_percentage")*d)/100.0);
    random_svd(A, A_rows_cols, A_rows_cols, A, L, V, A_rows_cols, k);

    return static_cast<unsigned long>(k);
}

template <typename T>
GurlsOptionsList *ParamSelHoPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //    [n,T]  = size(y);
    const unsigned long t = Y.cols();
    const unsigned long d = X.cols();


    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");


    const GurlsOptionsList* split = opt.getOptAs<GurlsOptionsList>("split");
    const gMat2D< unsigned long > &indices_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("indices");
    const gMat2D< unsigned long > &lasts_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("lasts");

    const unsigned long n = indices_mat.cols();

    const unsigned long *lasts = lasts_mat.getData();
    const unsigned long* indices_buffer = indices_mat.getData();


    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));

    int nholdouts = static_cast<int>(std::ceil( opt.getOptAsNumber("nholdouts")));

    gMat2D<T> *LAMBDA = new gMat2D<T>(1, t);
    T* lambdas = LAMBDA->getData();
    set(lambdas, (T)0.0, t);


    T* Q = new T[d*d];
    T* QtXty = new T[d*t];
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


    bool hasXt = opt.hasOpt("kernel.XtX") && opt.hasOpt("kernel.Xty");

    for(int nh=0; nh<nholdouts; ++nh)
    {
        unsigned long last = lasts[nh];
        unsigned long* tr = new unsigned long[last];
        unsigned long* va = new unsigned long[n-last];

        //copy int tr indices_ from n*nh to last
        copy< unsigned long >(tr,indices_buffer + n*nh,last,1,1);

        //copy int va indices_ from n*nh+last to n*nh+n
        copy< unsigned long >(va,(indices_buffer+ n*nh+last), n-last,1,1);


        gMat2D<T> Xva(n-last, d);
        gMat2D<T> yva(n-last, t);

        subMatrixFromRows(X.getData(), n, d, va, n-last, Xva.getData());
        subMatrixFromRows(Y.getData(), n, t, va, n-last, yva.getData());

        T* Xtr = NULL;
        if(hasXt)
        {
            T* XvatXva = new T[d*d];
            dot(Xva.getData(), Xva.getData(), XvatXva, n-last, d, n-last, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);

            const gMat2D<T>&XtX = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.XtX");

            // Q = XtX - XvatXva
            copy(Q, XtX.getData(), XtX.getSize());
            axpy(d*d, (T)-1.0, XvatXva, 1, Q, 1);

            delete [] XvatXva;
        }
        else
        {
            //       K = X(tr,:)'*X(tr,:);
            Xtr = new T[last*d];
            subMatrixFromRows(X.getData(), n, d, tr, last, Xtr);

            dot(Xtr, Xtr, Q, last, d, last, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
        }

        unsigned long k = eig_function(Q, L, d, d, opt, last);

        T* guesses = lambdaguesses(L, d, k, last, tot, (T)(opt.getOptAsNumber("smallnumber")));

        T* ap = new T[tot*t];


        if(hasXt)
        {
            T* Xvatyva = new T[d*t];
            dot(Xva.getData(), yva.getData(), Xvatyva, n-last, d, n-last, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

            const gMat2D<T>&Xty = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.Xty");


            // QtXty = Q'*(Xty - XvatXva)

            T* Xtrtytr = new T[d*t];

            copy(Xtrtytr, Xty.getData(), Xty.getSize());
            axpy(d*t, (T)-1.0, Xvatyva, 1, Xtrtytr, 1);

            dot(Q, Xtrtytr, QtXty, d, d, d, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

            delete [] Xvatyva;
            delete [] Xtrtytr;
        }
        else
        {
            T* ytr = new T[last*t];
            subMatrixFromRows(Y.getData(), n, t, tr, last, ytr);


            T* Xtrtytr = new T[d*t];
            dot(Xtr, ytr, Xtrtytr, last, d, last, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
            delete [] Xtr;

            dot(Q, Xtrtytr, QtXty, d, d, d, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

            delete [] ytr;
            delete [] Xtrtytr;
        }


        T* work = new T[d*(d+1)];

        for(int i=0; i<tot; ++i)
        {
            rls_eigen(Q, L, QtXty, W->getData(), guesses[i], last, d, d, d, d, t, work);

            OptMatrix<gMat2D<T> > *ret_pred = primal.execute(Xva, yva, *nestedOpt);

            nestedOpt->removeOpt("pred");
            nestedOpt->addOpt("pred", ret_pred);

            GurlsOptionsList* ret_perf = perfClass->execute(Xva, yva, *nestedOpt);

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
    delete [] QtXty;
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
