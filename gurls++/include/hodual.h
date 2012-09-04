

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
    virtual void eig_function(T* A, T* L, int A_rows_cols);
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
    virtual void eig_function(T* A, T* L, int A_rows_cols);
};





template<typename T>
void ParamSelHoDual<T>::eig_function(T* A, T* L, int A_rows_cols)
{
    eig_sm(A, L, A_rows_cols);
}

template<typename T>
void ParamSelHoDualr<T>::eig_function(T* A, T* L, int A_rows_cols)
{
    T* V = NULL;
    random_svd(A, A_rows_cols, A_rows_cols, A, L, V, A_rows_cols);
}

template <typename T>
GurlsOptionsList *ParamSelHoDual<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, const GurlsOptionsList &opt)
{
    //    [n,T]  = size(y);
    const int n = Y_OMR.rows();
    const int t = Y_OMR.cols();

    const int x_rows = X_OMR.rows();
    const int d = X_OMR.cols();

    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);


    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");
    nestedOpt->copyOpt<T>("kernel", opt);

    nestedOpt->addOpt("predkernel", new GurlsOptionsList("predkernel"));


    const GurlsOptionsList* split = GurlsOptionsList::dynacast(opt.getOpt("split"));
    const GurlsOption *index_opt = split->getOpt("indices");
    const GurlsOption *lasts_opt = split->getOpt("lasts");
    const gMat2D< unsigned long > &indices_mat = OptMatrix<gMat2D< unsigned long > >::dynacast(index_opt)->getValue();
    //vector of int one for column of indices_mat
    const gMat2D< unsigned long > &lasts_mat = OptMatrix<gMat2D< unsigned long > >::dynacast(lasts_opt)->getValue();
    const unsigned long *lasts_ = lasts_mat.getData();

    gMat2D< unsigned long > Indices(indices_mat.cols(), indices_mat.rows());
    indices_mat.transpose(Indices);
    unsigned long* indices_buffer = Indices.getData();

    const GurlsOptionsList* kernel = GurlsOptionsList::dynacast(opt.getOpt("kernel"));
    const gMat2D<T> &K = OptMatrix<gMat2D<T> >::dynacast(kernel->getOpt("K"))->getValue();

    const unsigned long k_rows = K.rows();

    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda") ));
    int nholdouts = static_cast<int>(std::ceil( opt.getOptAsNumber("nholdouts") ));

    T* lambdas = new T[t];
    set(lambdas, (T)0.0, t);

    T* acc_avg = new T[tot*t];
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
        unsigned long lasts = lasts_[nh];
        unsigned long* tr = new unsigned long[lasts];
        unsigned long* va = new unsigned long[n-lasts];

        //copy int tr indices_ from n*nh to lasts
        copy<unsigned long>(tr,indices_buffer + n*nh,lasts);
        //copy int va indices_ from n*nh+lasts to n*nh+n
        copy<unsigned long>(va,(indices_buffer+ n*nh+lasts), n-lasts);


        //Get K(tr,tr) from K
        T* Q = new T[lasts*lasts];
        copy_submatrix(Q, K.getData(), k_rows, lasts, lasts, tr,tr);

        T *L = new T[lasts];
        eig_function(Q, L, lasts);

        unsigned long r = lasts;

        if(kernel->getOptAsString("type") == "linear")
            r = std::min< unsigned long > (lasts,d);
        else
        {
            // 	opt.predkernel.K = opt.kernel.K(va,tr);%nva x ntr
            gMat2D<T> tmp(lasts, n-lasts);
            T* Kvatr = tmp.getData();
            copy_submatrix(Kvatr, K.getData(), k_rows, n-lasts, lasts, va, tr);

            GurlsOptionsList* predkernel = GurlsOptionsList::dynacast(nestedOpt->getOpt("predkernel"));

            gMat2D<T>* k_t = new gMat2D<T>(n-lasts, lasts);
            tmp.transpose(*k_t);
            predkernel->removeOpt("K");
            predkernel->addOpt("K", new OptMatrix<gMat2D<T> >(*k_t));
        }

        T* guesses = lambdaguesses(L, lasts, r, lasts, tot, (T)(opt.getOptAsNumber("smallnumber")));

        //  ap = zeros(tot,T);
        T* ap = new T[tot*t];
        T* Ytr = new T[lasts*t];
        subMatrixFromRows(Y.getData(), n, t, tr, lasts, Ytr);

        //    QtY = Q'*y(tr,:);
        T* Qty = new T[lasts*t];

        dot(Q, Ytr, Qty, lasts, lasts, lasts, t, lasts, t, CblasTrans, CblasNoTrans, CblasColMajor);

        T* C = new T[lasts*t];

        //    for i = 1:tot
        // 	opt.rls.X = X(tr,:);
        T* Xtr = new T[lasts*d];
        subMatrixFromRows(X.getData(), x_rows, d, tr, lasts, Xtr);

        gMat2D<T> rlsX_tmp(Xtr, d, lasts, false);
        gMat2D<T>* rlsX_t = new gMat2D<T>(lasts, d);
        rlsX_tmp.transpose(*rlsX_t);
        optimizer->removeOpt("X");
        optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*rlsX_t));

        T* Xva = new T[(n-lasts)*d];
        subMatrixFromRows(X.getData(), x_rows, d, va, n-lasts, Xva);

        T* Yva = new T[(n-lasts)*t];
        subMatrixFromRows(Y.getData(), n, t, va, n-lasts, Yva);

        for(int i=0; i<tot; ++i)
        {
            // 	opt.rls.C = rls_eigen(Q,L,QtY,guesses(i),ntr);
            rls_eigen(Q, L, Qty, C, guesses[i], lasts, lasts, lasts, lasts, lasts, t);

            gMat2D<T> c_tmp(C, t, lasts, false);
            gMat2D<T>* c_t = new gMat2D<T>(lasts, t);
            c_tmp.transpose(*c_t);

            optimizer->removeOpt("C");
            optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*c_t));

            if(kernel->getOptAsString("type") == "linear")
            {
//                opt.rls.W = X(tr,:)'*opt.rls.C; % dxT = (ntrxd)'*ntrxT   lasts*d lasts*lasts
                T* rlsW = new T[d*t];
                dot(Xtr, C, rlsW, lasts, d, lasts, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

                gMat2D<T> tmp(rlsW, t, d, false);
                gMat2D<T>* W = new gMat2D<T>(d, t);
                tmp.transpose(*W);

                delete[] rlsW;
                optimizer->removeOpt("W");
                optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

            }
//            else
//            opt.predkernel.K = opt.kernel.K(va,tr);%nva x ntr
//            (see above)


            gMat2D<T> tmp_x(Xva, d, n -lasts, false);
            gMat2D<T>* xx = new gMat2D<T>(n-lasts, d);
            tmp_x.transpose(*xx);


            gMat2D<T> tmp_y(Yva,t, n-lasts, false);
            gMat2D<T>* yy = new gMat2D<T>(n-lasts, t);
            tmp_y.transpose(*yy);

            OptMatrix<gMat2D<T> >* ret_pred = dual.execute(*xx, *yy, *nestedOpt);

            nestedOpt->addOpt("pred", ret_pred);

            // 	opt.perf = opt.hoperf(Xva,yva,opt);
            GurlsOptionsList* ret_perf = perfClass->execute(*xx, *yy, *nestedOpt);

            gMat2D<T> &forho_vec = OptMatrix<gMat2D<T> >::dynacast(ret_perf->getOpt("forho"))->getValue();

            //       for t = 1:T
            //          ap(i,t) = opt.perf.forho(t);
            copy(ap+i, forho_vec.getData(), t, tot, 1);

            nestedOpt->removeOpt("pred");
            delete ret_perf;

        }//for tot

        delete [] Q;
        delete [] Qty;
        delete [] L;
        delete [] C;
        delete [] va;
        delete [] tr;
        delete [] Ytr;
        delete [] Yva;
        delete [] Xva;
        delete [] Xtr;

        //[dummy,idx] = max(ap,[],1);
        T* work = NULL;
        unsigned long* idx = new unsigned long[t];
        indicesOfMax(ap, tot, t, idx, work, 1);

//        delete[] work;


        //vout.lambdas_round{nh} = guesses(idx);
        T* lambdas_nh = copyLocations(idx, guesses, t, tot);

        copy(lambdas_round+nh*t, lambdas_nh, t);

        //add lambdas_nh to lambdas
        axpy< T >(t, (T)1, lambdas_nh, 1, lambdas, 1);

        delete [] lambdas_nh;
        delete [] idx;

        //  vout.perf{nh} = ap;
        copy(perf + nh*tot*t, ap, tot*t, 1, 1);

        axpy(tot*t, (T)1, ap, 1, acc_avg, 1);

        //  vout.guesses{nh} = guesses;
        copy(ret_guesses + nh*tot,guesses,tot,1,1);

        delete [] guesses;
        delete [] ap;
    }//for nholdouts

    delete nestedOpt;
    delete perfClass;


    GurlsOptionsList* paramsel;

    if(opt.hasOpt("paramsel"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt<T>("paramsel", opt);

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

    gMat2D< T > acc_tmp(acc_avg,t,tot,false);
    gMat2D<T>*  acc_t = new gMat2D<T>(tot, t);
    acc_tmp.transpose(*acc_t);
    delete [] acc_avg;
    paramsel->addOpt("acc_avg", new OptMatrix<gMat2D<T> >(*acc_t));


    gMat2D<T> *LAMBDA = new gMat2D<T>(1, t);
    copy(LAMBDA->getData(), lambdas, t);

    delete [] lambdas;

    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*LAMBDA));
    paramsel->addOpt("lambdas_round", new OptMatrix<gMat2D<T> >(*lambdas_round_mat));


    return paramsel;

}


}

#endif // _GURLS_HODUAL_H_
