

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


#ifndef _GURLS_HODUAL_H_
#define _GURLS_HODUAL_H_

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

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelHoDual is the subclass of ParamSelection that implements hold-out cross validation with the dual formulation of RLS
 */

template <typename T>
class ParamSelHoDual: public ParamSelection<T>{

public:
    /**
     * Performs parameter selection when the dual formulation of RLS is used.
     * The hold-out approach is used.
     * The performance measure specified by opt.hoperf is maximized.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - nlambda (default)
     *  - hoperf (default)
     *  - smallnumber (default)
     *  - split (settable with the class Split and its subclasses)
     *  - kernel (settable with the class Kernel and its subclasses)
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
    virtual void eig_function(T* A, T* L, int A_rows_cols, unsigned long n, const GurlsOptionsList &opt);
};


/**
 * \ingroup ParameterSelection
 * \brief ParamSelHoDualr is the randomized version of \ref ParamSelHoDual
 */

template <typename T>
class ParamSelHoDualr: public ParamSelHoDual<T>{

protected:
    /**
     * Auxiliary method used to call the right eig/svd function for this class
     */
    virtual void eig_function(T* A, T* L, int A_rows_cols, unsigned long n, const GurlsOptionsList &opt);
};





template<typename T>
void ParamSelHoDual<T>::eig_function(T* A, T* L, int A_rows_cols, unsigned long , const GurlsOptionsList &)
{
    eig_sm(A, L, A_rows_cols);
}

template<typename T>
void ParamSelHoDualr<T>::eig_function(T* A, T* L, int A_rows_cols, unsigned long n, const GurlsOptionsList &opt)
{
    T* V = NULL;
    T k = gurls::round((opt.getOptAsNumber("eig_percentage")*n)/100.0);
    random_svd(A, A_rows_cols, A_rows_cols, A, L, V, k);
}

template <typename T>
GurlsOptionsList *ParamSelHoDual<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //    [n,T]  = size(y);
    const unsigned long n = Y.rows();
    const unsigned long t = Y.cols();

    const unsigned long x_rows = X.rows();
    const unsigned long d = X.cols();


    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");
    nestedOpt->copyOpt("kernel", opt);

    nestedOpt->addOpt("predkernel", new GurlsOptionsList("predkernel"));


    const GurlsOptionsList* split = opt.getOptAs<GurlsOptionsList>("split");

    const gMat2D< unsigned long > &indices_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("indices");
    const gMat2D< unsigned long > &lasts_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("lasts");

    const unsigned long *lasts = lasts_mat.getData();
    const unsigned long* indices_buffer = indices_mat.getData();


    const GurlsOptionsList* kernel = opt.getOptAs<GurlsOptionsList>("kernel");
    const gMat2D<T> &K = kernel->getOptValue<OptMatrix<gMat2D<T> > >("K");
    const bool linearKernel = kernel->getOptAsString("type") == "linear";

    const unsigned long k_rows = K.rows();

    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));
    int nholdouts = static_cast<int>(std::ceil( opt.getOptAsNumber("nholdouts")));

    gMat2D<T> *LAMBDA = new gMat2D<T>(1, t);
    T* lambdas = LAMBDA->getData();
    set(lambdas, (T)0.0, t);

    gMat2D<T>*  acc_avg_mat = new gMat2D<T>(tot, t);
    T* acc_avg = acc_avg_mat->getData();
    set(acc_avg, (T)0.0, tot*t);

    //     for nh = 1:opt.nholdouts
    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));
    PredDual<T> dual;

    nestedOpt->addOpt("optimizer",optimizer);

    gMat2D<T>* perf_mat = new gMat2D<T>(nholdouts, tot*t);
    T* perf = perf_mat->getData();

    gMat2D<T>*  guesses_mat = new gMat2D<T>(nholdouts, tot);
    T* ret_guesses = guesses_mat->getData();

    gMat2D<T>* lambdas_round_mat = new gMat2D<T>(nholdouts, t);
    T* lambdas_round = lambdas_round_mat->getData();

    for(int nh=0; nh<nholdouts; ++nh)
    {
        unsigned long last = lasts[nh];
        unsigned long* tr = new unsigned long[last];
        unsigned long* va = new unsigned long[n-last];

        //copy int tr indices_ from n*nh to last
        copy<unsigned long>(tr,indices_buffer + n*nh, last);
        //copy int va indices_ from n*nh+last to n*nh+n
        copy<unsigned long>(va,(indices_buffer+ n*nh+last), n-last);


        //Get K(tr,tr) from K
        T* Q = new T[last*last];
        copy_submatrix(Q, K.getData(), k_rows, last, last, tr, tr);

        T *L = new T[last];
        eig_function(Q, L, last, n, opt);

        unsigned long r = last;

        if(linearKernel)
            r = std::min< unsigned long > (last,d);
        else
        {
            // 	opt.predkernel.K = opt.kernel.K(va,tr);%nva x ntr
            gMat2D<T>* predK = new gMat2D<T>(n-last, last);

            copy_submatrix(predK->getData(), K.getData(), k_rows, n-last, last, va, tr);

            GurlsOptionsList* predkernel = nestedOpt->getOptAs<GurlsOptionsList>("predkernel");

            predkernel->removeOpt("K");
            predkernel->addOpt("K", new OptMatrix<gMat2D<T> >(*predK));
        }

        T* guesses = lambdaguesses(L, last, r, last, tot, (T)(opt.getOptAsNumber("smallnumber")));

        //  ap = zeros(tot,T);
        T* ap = new T[tot*t];

        T* Ytr = new T[last*t];
        subMatrixFromRows(Y.getData(), n, t, tr, last, Ytr);

        //    QtY = Q'*y(tr,:);
        T* Qty = new T[last*t];

        dot(Q, Ytr, Qty, last, last, last, t, last, t, CblasTrans, CblasNoTrans, CblasColMajor);

        delete [] Ytr;

        //    for i = 1:tot
        // 	opt.rls.X = X(tr,:);

        gMat2D<T>* rlsX = new gMat2D<T>(last, d);
        subMatrixFromRows(X.getData(), x_rows, d, tr, last, rlsX->getData());

        delete [] tr;

        optimizer->removeOpt("X");
        optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*rlsX));


        gMat2D<T>* xx = new gMat2D<T>(n-last, d);
        subMatrixFromRows(X.getData(), x_rows, d, va, n-last, xx->getData());

        gMat2D<T>* yy = new gMat2D<T>(n-last, t);
        subMatrixFromRows(Y.getData(), n, t, va, n-last, yy->getData());

        delete [] va;

        gMat2D<T>* C = new gMat2D<T>(last, t);
        optimizer->removeOpt("C");
        optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*C));


        gMat2D<T>* W = NULL;

        if(linearKernel)
        {
            W = new gMat2D<T>(d, t);

            optimizer->removeOpt("W");
            optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));
        }

        T* work = new T[last*(last+1)];

        for(int i=0; i<tot; ++i)
        {
            // 	opt.rls.C = rls_eigen(Q,L,QtY,guesses(i),ntr);
            rls_eigen(Q, L, Qty, C->getData(), guesses[i], last, last, last, last, last, t, work);

            if(linearKernel)
            {
//                opt.rls.W = X(tr,:)'*opt.rls.C; % dxT = (ntrxd)'*ntrxT   last*d last*last
                dot(rlsX->getData(), C->getData(), W->getData(), last, d, last, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
            }
//            else
//            opt.predkernel.K = opt.kernel.K(va,tr);%nva x ntr
//            (see above)


            OptMatrix<gMat2D<T> >* ret_pred = dual.execute(*xx, *yy, *nestedOpt);

            nestedOpt->addOpt("pred", ret_pred);

            // 	opt.perf = opt.hoperf(Xva,yva,opt);
            GurlsOptionsList* ret_perf = perfClass->execute(*xx, *yy, *nestedOpt);

            gMat2D<T> &forho_vec = ret_perf->getOptValue<OptMatrix<gMat2D<T> > >("forho");

            //       for t = 1:T
            //          ap(i,t) = opt.perf.forho(t);
            copy(ap+i, forho_vec.getData(), t, tot, 1);

            nestedOpt->removeOpt("pred");

            delete ret_perf;

        }//for tot

        delete [] Q;
        delete [] Qty;
        delete [] L;
        delete [] work;

        delete xx;
        delete yy;

        //[dummy,idx] = max(ap,[],1);
        work = NULL;
        unsigned long* idx = new unsigned long[t];
        indicesOfMax(ap, tot, t, idx, work, 1);


        //vout.lambdas_round{nh} = guesses(idx);
        T* lambdas_nh = new T[t];
        copyLocations(idx, guesses, t, tot, lambdas_nh);

        copy(lambdas_round +nh, lambdas_nh, t, nholdouts, 1);

        //add lambdas_nh to lambdas
        axpy< T >(t, (T)1, lambdas_nh, 1, lambdas, 1);

        delete [] lambdas_nh;
        delete [] idx;

        //  vout.perf{nh} = ap;
        copy(perf + nh, ap, tot*t, nholdouts, 1);

        axpy(tot*t, (T)1, ap, 1, acc_avg, 1);

        //  vout.guesses{nh} = guesses;
        copy(ret_guesses + nh, guesses, tot, nholdouts, 1);

        delete [] guesses;
        delete [] ap;
    }//for nholdouts

    delete nestedOpt;
    delete perfClass;


    GurlsOptionsList* paramsel;

    if(opt.hasOpt("paramsel"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt("paramsel", opt);

        paramsel = GurlsOptionsList::dynacast(tmp_opt->getOpt("paramsel"));
        tmp_opt->removeOpt("paramsel", false);
        delete tmp_opt;

        paramsel->removeOpt("perf");
        paramsel->removeOpt("guesses");
        paramsel->removeOpt("acc_avg");
        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("lambdas_round");
    }
    else
        paramsel = new GurlsOptionsList("paramsel");



    paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*perf_mat));
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

    //     if numel(vout.lambdas_round) > 1
    // 	lambdas = cell2mat(vout.lambdas_round');
    // 	vout.lambdas = mean(lambdas);
    //     else
    // 	vout.lambdas = vout.lambdas_round{1};
    //     end
    if(nholdouts>1)
    {
        scal(t, (T)1.0/nholdouts, lambdas, 1);
        scal(t, (T)1.0/nholdouts, acc_avg, 1);
    }

    paramsel->addOpt("acc_avg", new OptMatrix<gMat2D<T> >(*acc_avg_mat));


    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*LAMBDA));
    paramsel->addOpt("lambdas_round", new OptMatrix<gMat2D<T> >(*lambdas_round_mat));

    return paramsel;
}


}

#endif // _GURLS_HODUAL_H_
