

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
   * \brief HoDual is the subclass of ParamSelection that implements hold-out cross validation with the dual formulation of RLS
   */

template <typename T>
class HoDual: public ParamSelection<T>{

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
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void HoDual<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
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

    GurlsOptionsList* predkernel_old = NULL;
    if(opt.hasOpt("predkernel"))
    {
        predkernel_old = static_cast<GurlsOptionsList*>(opt.getOpt("predkernel"));
        opt.removeOpt("predkernel", false);
    }
    opt.addOpt("predkernel", new GurlsOptionsList("predkernel"));


    GurlsOptionsList* split = static_cast<GurlsOptionsList*>(opt.getOpt("split"));
    GurlsOption *index_opt = split->getOpt("indices");
    GurlsOption *lasts_opt = split->getOpt("lasts");
    gMat2D< unsigned long > *indices_mat = &(OptMatrix<gMat2D< unsigned long > >::dynacast(index_opt))->getValue();
    //vector of int one for column of indices_mat
    gMat2D< unsigned long > *lasts_mat = &(OptMatrix<gMat2D< unsigned long > >::dynacast(lasts_opt))->getValue();
    unsigned long *lasts_ = lasts_mat->getData();

    gMat2D< unsigned long > Indices(indices_mat->cols(), indices_mat->rows());
    indices_mat->transpose(Indices);
    unsigned long* indices_buffer = new unsigned long[indices_mat->cols()*indices_mat->rows()];
    copy< unsigned long >(indices_buffer,Indices.getData(),indices_mat->cols()*indices_mat->rows());

    GurlsOptionsList* kernel = static_cast<GurlsOptionsList*>(opt.getOpt("kernel"));
    GurlsOption *K_opt = kernel->getOpt("K");
    gMat2D<T> *K_mat = &(OptMatrix<gMat2D<T> >::dynacast(K_opt))->getValue();
    gMat2D<T> K(K_mat->cols(), K_mat->rows());
    K_mat->transpose(K);
    unsigned long k_rows = K_mat->rows();
    T* k_buffer_total = new T[K_mat->cols()*K_mat->rows()];
    copy< T >(k_buffer_total,K.getData(),K_mat->cols()*K_mat->rows());

    unsigned long* indices_from0toD = new unsigned long[d];
    unsigned long* indices_from0toT = new unsigned long[t];

    for(int i=0;i<d;++i)
        indices_from0toD[i]=i;
    for(int i=0;i<t;++i)
        indices_from0toT[i]=i;

    int tot = static_cast<int>(std::ceil( static_cast<OptNumber*>(opt.getOpt("nlambda"))->getValue() ));
    int nholdouts = static_cast<int>(std::ceil( static_cast<OptNumber*>(opt.getOpt("nholdouts"))->getValue() ));

    T* lambdas = new T[t];
    set(lambdas, (T)0.0, t);

    T* forho_nholdouts = new T[nholdouts * tot*t];
    T* acc_avg = new T[tot*t];
    set(acc_avg, (T)0.0, tot*t);
    T* guesses_nholdouts = new T[nholdouts * tot];

    //     for nh = 1:opt.nholdouts
    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");
    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));
    PredDual< T > dual;

    opt.addOpt("optimizer",optimizer);

    for(int nh=0; nh<nholdouts; ++nh)
    {
        unsigned long lasts = lasts_[nh];
        unsigned long* tr = new unsigned long[lasts];
        unsigned long* va = new unsigned long[n-lasts];

        //copy int tr indices_ from n*nh to lasts
        copy<unsigned long>(tr,indices_buffer + n*nh,lasts);
        //copy int va indices_ from n*nh+lasts to n*nh+n
        copy<unsigned long>(va,(indices_buffer+ n*nh+lasts), n-lasts);


        //Get K(tr,tr) from k_buffer_total
        T* Q = new T[lasts*lasts];
        copy_submatrix(Q, k_buffer_total, k_rows, lasts, lasts, tr,tr);

        T *L = new T[lasts];
        eig_sm(Q, L, lasts, lasts);

        unsigned long r = lasts;
        if(kernel->getOptAsString("type") == "linear")
            r = std::min< unsigned long > (lasts,d);
        else
        {
            // 	opt.predkernel.K = opt.kernel.K(va,tr);%nva x ntr
            gMat2D<T> tmp(lasts, n-lasts);
            T* Kvatr = tmp.getData();
            copy_submatrix(Kvatr, k_buffer_total, k_rows, n-lasts, lasts, va, tr);

            GurlsOptionsList* predkernel = static_cast<GurlsOptionsList*>(opt.getOpt("predkernel"));

            gMat2D<T>* k_t = new gMat2D<T>(n-lasts, lasts);
            tmp.transpose(*k_t);
            predkernel->removeOpt("K");
            predkernel->addOpt("K", new OptMatrix<gMat2D<T> >(*k_t));
        }

        T* guesses = lambdaguesses(L, lasts, r, lasts, tot, (T)(opt.getOptAsNumber("smallnumber")));

        //  ap = zeros(tot,T);
        T* ap = new T[tot*t];
        T* Ytr = new T[lasts*t];
        copy_submatrix(Ytr, Y.getData(), n, lasts, t, tr, indices_from0toT);
        //    QtY = Q'*y(tr,:); %
        T* Qty = new T[lasts*t];

        //       dot(Q, Y.getData(), Qty, lasts, lasts, lasts, t, lasts, t, CblasTrans, CblasNoTrans, CblasColMajor);
        dot(Q, Ytr, Qty, lasts, lasts, lasts, t, lasts, t, CblasTrans, CblasNoTrans, CblasColMajor);

        T* C = new T[lasts*t];
        //    for i = 1:tot
        // 	opt.rls.X = X(tr,:);
        T* Xtr = new T[lasts*d];
        copy_submatrix(Xtr, X.getData(), x_rows, lasts, d, tr, indices_from0toD);
        gMat2D<T> rlsX_tmp(Xtr, d, lasts, false);
        gMat2D<T>* rlsX_t = new gMat2D<T>(lasts, d);
        rlsX_tmp.transpose(*rlsX_t);
        optimizer->removeOpt("X");
        optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*rlsX_t));

        T* Xva = new T[(n-lasts)*d];
        copy_submatrix(Xva, X.getData(), x_rows, n-lasts, d, va, indices_from0toD);

        T* Yva = new T[(n-lasts)*t];
        copy_submatrix(Yva, Y.getData(), n, n-lasts, t, va, indices_from0toT);

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
                //         opt.rls.W = X(tr,:)'*opt.rls.C; % dxT = (ntrxd)'*ntrxT   lasts*d lasts*lasts
                T* rlsW = new T[d*t];
                dot(Xtr, C, rlsW, lasts, d, lasts, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
                gMat2D<T> tmp(rlsW, t, d, false);
                gMat2D<T>* W = new gMat2D<T>(d, t);
                tmp.transpose(*W);
                delete[] rlsW;
                optimizer->removeOpt("W");
                optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

            }

            // 	%elseif strcmp(opt.kernel.type,'load')
            // 	%	opt.predkernel.type = 'load';
            // 	else
//            else
//            {
//                // 	opt.predkernel.K = opt.kernel.K(va,tr);%nva x ntr
//                T* Kvatr = new T[(n-lasts)*lasts];
//                copy_submatrix(Kvatr, k_buffer_total, k_rows, n-lasts, lasts, va, tr);
//                if(opt.hasOpt())
//                GurlsOptionsList* predkernel = static_cast<GurlsOptionsList*>(opt.getOpt("predkernel"));

//                gMat2D<T> tmp(Kvatr, lasts, n-lasts, false);
//                gMat2D<T>* k_t = new gMat2D<T>(n-lasts, lasts);
//                tmp.transpose(*k_t);
//                optimizer->removeOpt("K");
//                predkernel->addOpt("K", new OptMatrix<gMat2D<T> >(*k_t));
//                delete[] Kvatr;
//            }
            // 	%else
            // 	%	opt.predkernel = predkernel_traintest(X(va,:),y(va,:),opt);
            //      end
            // 	opt.pred = pred_dual(Xva,yva,opt);


            gMat2D<T> tmp_x(Xva, d, n -lasts, false);
            gMat2D<T>* xx = new gMat2D<T>(n-lasts, d);
            tmp_x.transpose(*xx);


            gMat2D<T> tmp_y(Yva,t, n-lasts, false);
            gMat2D<T>* yy = new gMat2D<T>(n-lasts, t);
            tmp_y.transpose(*yy);

            dual.execute(*xx, *yy, opt);


            //  	gMat2D<T> *predV = &(OptMatrix<gMat2D<T> >::dynacast(opt.getOpt("pred")))->getValue();

            // 	opt.perf = opt.hoperf(Xva,yva,opt);
            perfClass->execute(*xx, *yy, opt);

            GurlsOptionsList* perf = static_cast<GurlsOptionsList*>(opt.getOpt("perf"));

            gMat2D<T> *forho_vec = &(OptMatrix<gMat2D<T> >::dynacast(perf->getOpt("forho")))->getValue();
            //       for t = 1:T
            for(int j = 0; j<t; ++j)
            {
                //          ap(i,t) = opt.perf.forho(t);
                ap[i +(tot*j)] = forho_vec->getData()[j];
            }

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

        delete[] work;


        //vout.lambdas_round{nh} = guesses(idx);
        T* lambdas_nh = copyLocations(idx, guesses, t, tot);
        //add lambdas_nh to lambdas
        axpy< T >(t, (T)1, lambdas_nh, 1, lambdas, 1);
        delete [] lambdas_nh;
        delete [] idx;
        //  vout.forho{nh} = ap;
        copy< T >(forho_nholdouts + nh*tot*t,ap,tot*t,1,1);
    axpy< T >(tot*t, (T)1, ap, 1, acc_avg, 1);
        //  vout.guesses{nh} = guesses;
        copy< T >(guesses_nholdouts + nh*tot,guesses,tot,1,1);

        delete [] guesses;
        delete [] ap;
    }//for nholdouts

    delete [] indices_from0toD;
    delete [] indices_from0toT;

    delete perfClass;

    delete [] indices_buffer;
    delete [] k_buffer_total;

    gMat2D< T > forho_tmp(forho_nholdouts, tot*t,nholdouts,false);
    gMat2D<T>* forho_t = new gMat2D<T>(nholdouts, tot*t);
    forho_tmp.transpose(*forho_t);

    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");
    paramsel->addOpt("forho", new OptMatrix<gMat2D<T> >(*forho_t));

    delete [] forho_nholdouts;

    gMat2D< T > guesses_tmp(guesses_nholdouts,tot,nholdouts,false);
    gMat2D<T>*  guesses_t = new gMat2D<T>(nholdouts, tot);
    guesses_tmp.transpose(*guesses_t);
    //     opt.addOpt("guesses", new OptMatrix<gMat2D<T> >(guesses_t));
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_t));

    delete [] guesses_nholdouts;

    //     if numel(vout.lambdas_round) > 1
    // 	lambdas = cell2mat(vout.lambdas_round');
    // 	vout.lambdas = mean(lambdas);
    //     else
    // 	vout.lambdas = vout.lambdas_round{1};
    //     end
    if(nholdouts>1)
    {
        for(T *it = lambdas, *end = lambdas+t; it != end; ++it)
            *it = *it / nholdouts;
    for(T *it = acc_avg, *end = acc_avg+tot*t; it != end; ++it)
            *it = *it / nholdouts;
    }
    gMat2D< T > acc_tmp(acc_avg,t,tot,false);
    gMat2D<T>*  acc_t = new gMat2D<T>(tot, t);
    acc_tmp.transpose(*acc_t);
    paramsel->addOpt("acc_avg", new OptMatrix<gMat2D<T> >(*acc_t));

    OptNumberList* LAMBDA = new OptNumberList();
    for (T* l_it = lambdas, *l_end = lambdas+t; l_it != l_end; ++l_it)
    {
        LAMBDA->add(static_cast<double>(*l_it));
    }

    delete [] lambdas;

    //     opt.addOpt("lambdas", LAMBDA);
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

    opt.removeOpt("predkernel");
    if(predkernel_old != NULL)
        opt.addOpt("predkernel", predkernel_old);

}


}

#endif // _GURLS_HODUAL_H_
