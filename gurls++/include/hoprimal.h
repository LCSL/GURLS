

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

namespace gurls {

/**
   * \brief HoPrimal is the subclass of ParamSelection that implements hold-out cross validation with the primal formulation of RLS
   */

template <typename T>
class HoPrimal: public ParamSelection<T>{

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
};

template <typename T>
void HoPrimal<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
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
        perf_old = static_cast<GurlsOptionsList*>(opt.getOpt("perf"));
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
        optimizer_old = static_cast<GurlsOptionsList*>(opt.getOpt("optimizer"));
        opt.removeOpt("optimizer", false);
    }

    GurlsOptionsList* split = static_cast<GurlsOptionsList*>(opt.getOpt("split"));
    GurlsOption *index_opt = split->getOpt("indices");
    GurlsOption *lasts_opt = split->getOpt("lasts");

    GurlsOption *index_monline_opt;
    gMat2D< unsigned long > *indices_monline_mat;
    unsigned long* indices_monline_buffer = NULL;
    if(opt.hasOpt("hoMOnline"))
    {
        //tr = opt.split.monline{1}.idx;
        index_monline_opt = split->getOpt("monline");
        indices_monline_mat = &(OptMatrix<gMat2D< unsigned long > >::dynacast(index_monline_opt))->getValue();
        gMat2D< unsigned long > Indices(indices_monline_mat->cols(), indices_monline_mat->rows());
        indices_monline_mat->transpose(Indices);
        indices_monline_buffer = new unsigned long[indices_monline_mat->rows()];
        copy< unsigned long >(indices_monline_buffer,Indices.getData(),indices_monline_mat->rows());
    }

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

    int tot = static_cast<int>(std::ceil( static_cast<OptNumber*>(opt.getOpt("nlambda"))->getValue() ));

    int nholdouts = static_cast<int>(std::ceil( static_cast<OptNumber*>(opt.getOpt("nholdouts"))->getValue() ));
    T* lambdas = new T[t];
    set(lambdas, (T)0.0, t);

    T* forho_nholdouts = new T[nholdouts * tot*t];
    T* guesses_nholdouts = new T[nholdouts * tot];

    T* Q = new T[d*d];
    T* Qty = new T[d*t];
    T *L = new T[d];
    T* W = new T[d*t];

    for(int nh=0;nh<nholdouts;nh++)
    {
        unsigned long lasts = lasts_[nh];
        unsigned long* tr = new unsigned long[lasts];
        unsigned long* va = new unsigned long[n-lasts];

        if(opt.hasOpt("hoMOnline"))
            copy< unsigned long >(tr,indices_monline_buffer,lasts,1,1);
        else
            //copy int tr indices_ from n*nh to lasts
            copy< unsigned long >(tr,indices_buffer + n*nh,lasts,1,1);

        //copy int va indices_ from n*nh+lasts to n*nh+n
        copy< unsigned long >(va,(indices_buffer+ n*nh+lasts), n-lasts,1,1);

        const int qrows = d;
        const int qcols = d;
        //       const int l_length = qrows;
        //       K = X(tr,:)'*X(tr,:);
        T* Xtr = new T[lasts*d];
        copy_submatrix(Xtr, X.getData(), x_rows, lasts, d, tr, indices_from0toD);
        //Get K(tr,tr) from k_buffer_total
        //       T* Q = new T[d,d];
        dot(Xtr, Xtr, Q, lasts, d, lasts, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
        //       T *L = new T[l_length];
        eig_sm(Q, L, qrows, qcols);

        unsigned long miN = d<lasts?d:lasts;
        T* guesses = lambdaguesses(L, d, miN,lasts, tot, (T)(opt.getOptAsNumber("smallnumber")));

        T* ap = new T[tot*t];

        T* Ytr = new T[lasts*t];
        copy_submatrix(Ytr, Y.getData(), n, lasts, t, tr, indices_from0toT);


        T* QtyT = new T[d * lasts];
        dot(Q, Xtr, QtyT, d, d, lasts, d, d, lasts, CblasTrans, CblasTrans, CblasColMajor);

        dot(QtyT, Ytr, Qty, qrows, lasts, lasts, t, qcols, t, CblasNoTrans, CblasNoTrans, CblasColMajor);
        delete [] QtyT;

        GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");
        opt.addOpt("optimizer",optimizer);

        for(int i=0;i<tot;i++)
        {
            rls_eigen(Q, L, Qty, W, guesses[i], lasts, d, d, d, qcols, t);

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
            gMat2D<T>* xx = new gMat2D<T>(n-lasts, d);
            Xva_m.transpose(*xx);


            const gMat2D<T> Yva_m(Yva, t, n-lasts, false);
            gMat2D<T>* yy = new gMat2D<T>(n-lasts, t);
            Yva_m.transpose(*yy);


            PredPrimal< T > primal;

            primal.execute(*xx, *yy, opt);

            Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));
            perfClass->execute( *xx, *yy, opt);

            delete perfClass;
            delete [] Yva;
            delete [] Xva;

            GurlsOptionsList* perf = static_cast<GurlsOptionsList*>(opt.getOpt("perf"));
            gMat2D<T> *forho_vec = &(OptMatrix<gMat2D<T> >::dynacast(perf->getOpt("forho")))->getValue();
            for(unsigned long j = 0; j<t; ++j)
            {
                ap[i +(tot*j)] = forho_vec->getData()[j];
            }
        }

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
        //add lambdas_nh to lambdas
        axpy< T >(t, (T)1, lambdas_nh, 1, lambdas, 1);
        delete [] lambdas_nh;
        delete [] idx;

        //  vout.forho{nh} = ap;
        copy< T >(forho_nholdouts + nh*tot*t,ap,tot*t,1,1);
        //  vout.guesses{nh} = guesses;
        copy< T >(guesses_nholdouts + nh*tot,guesses,tot,1,1);
        
        delete [] guesses;
        delete [] ap;

    }

    if(indices_monline_buffer!=NULL)
        delete [] indices_monline_buffer;

    delete [] indices_buffer;
    delete [] Q;
    delete [] Qty;
    delete [] L;
    delete [] W;


    gMat2D< T > forho_tmp(forho_nholdouts,tot*t,nholdouts,false);
    gMat2D<T>* forho_t = new gMat2D<T>(nholdouts, tot*t);
    forho_tmp.transpose(*forho_t);
    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");
    paramsel->addOpt("forho", new OptMatrix<gMat2D<T> >(*forho_t));

    delete [] forho_nholdouts;

    gMat2D< T > guesses_tmp(guesses_nholdouts,tot,nholdouts,false);
    gMat2D<T>*  guesses_t = new gMat2D<T>(nholdouts, tot);
    guesses_tmp.transpose(*guesses_t);
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_t));

    delete [] guesses_nholdouts;

    if(nholdouts>1)
        for(T *it = lambdas, *end = lambdas+t; it != end; ++it)
            *it = *it / nholdouts;

    OptNumberList* LAMBDA = new OptNumberList();
    for (T* l_it = lambdas, *l_end = lambdas+t; l_it != l_end; ++l_it)
        LAMBDA->add(static_cast<double>(*l_it));

    delete [] lambdas;

    paramsel->addOpt("lambdas", LAMBDA);
    opt.addOpt("paramsel", paramsel);

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
