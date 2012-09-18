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


#ifndef _GURLS_HOGPREGR_H_
#define _GURLS_HOGPREGR_H_

#include <cmath>

#include "options.h"
#include "optlist.h"
#include "gmat2d.h"
#include "gvec.h"
#include "gmath.h"

#include "paramsel.h"
#include "perf.h"

#include "rlsgp.h"
#include "predgp.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelHoGPRegr is the sub-class of ParamSelection that implements
 */

template <typename T>
class ParamSelHoGPRegr: public ParamSelection<T>{

public:
    /**
     *
     */
   GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList *ParamSelHoGPRegr<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //    [n,T]  = size(y);
    const unsigned long n = Y.rows();
    const unsigned long t = Y.cols();

    const unsigned long d = X.cols();

//    tot = opt.nlambda;
    int tot = static_cast<int>(opt.getOptAsNumber("nlambda"));


//    K = opt.kernel.K;
    const gMat2D<T> &K = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.K");

    const GurlsOptionsList* split = opt.getOptAs<GurlsOptionsList>("split");

    const gMat2D< unsigned long > &indices_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("indices");
    const gMat2D< unsigned long > &lasts_mat = split->getOptValue<OptMatrix<gMat2D< unsigned long > > >("lasts");

    const unsigned long *lasts = lasts_mat.getData();
    const unsigned long *indices = indices_mat.getData();


    unsigned long *tr = new unsigned long[indices_mat.cols()];
    unsigned long *va;


    T* work = new T[t+n];

//    lmax = mean(std(y));
    T* stdY = new T[t];

    stdDev(Y.getData(), n, t, stdY, work);

    const T lmax = sumv(stdY, t, work)/((T)t);

//    lmin = mean(std(y))*10^-5;
    const T lmin = lmax * (T)1.0e-5;

    delete[] stdY;
    delete[] work;

//    guesses = lmin.*(lmax/lmin).^linspace(0,1,tot);
    T* guesses = new T[tot];

    T* linspc = new T[tot];
    linspace((T)0.0, (T)1.0, tot, linspc);
    const T coeff = lmax/lmin;

    for(int i=0; i< tot; ++i)
        guesses[i] = lmin* std::pow(coeff, linspc[i]);

    delete[] linspc;

    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");
//    nestedOpt->copyOpt<T>("singlelambda", opt);


    GurlsOptionsList* tmpPredKernel = new GurlsOptionsList("predkernel");
    GurlsOptionsList* tmpKernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* tmpParamSel = new GurlsOptionsList("paramsel");

    nestedOpt->addOpt("kernel", tmpKernel);
    nestedOpt->addOpt("predkernel", tmpPredKernel);
    nestedOpt->addOpt("paramsel", tmpParamSel);

    gMat2D<T> subXtr;
    gMat2D<T> subYtr;
    gMat2D<T> subXva;
    gMat2D<T> subYva;

    gMat2D<T>* subK = new gMat2D<T>();
    gMat2D<T>* subPredK = new gMat2D<T>();
    gMat2D<T>* subPredKTest = new gMat2D<T>();

    tmpKernel->addOpt("K", new OptMatrix<gMat2D<T> > (*subK));
    tmpPredKernel->addOpt("K", new OptMatrix<gMat2D<T> > (*subPredK));
    tmpPredKernel->addOpt("Ktest", new OptMatrix<gMat2D<T> > (*subPredKTest));


    RLSGPRegr<T> rlsgp;
    PredGPRegr<T> predgp;
    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));

    const int nholdouts = static_cast<int>(opt.getOptAsNumber("nholdouts"));
    const unsigned long indices_cols = indices_mat.cols();

    T *perf = new T[tot*t];

    gMat2D<T>* lambdas_round_mat = new gMat2D<T>(nholdouts, t);
    T *lambdas_round = lambdas_round_mat->getData();

    gMat2D<T> *perf_mat = new gMat2D<T>(nholdouts, tot*t);

    gMat2D<T>* guesses_mat = new gMat2D<T>(nholdouts, tot);
    T *ret_guesses = guesses_mat->getData();

    gMat2D<T> * lambda = new gMat2D<T>(1,1);
    tmpParamSel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*lambda));

//    for nh = 1:opt.nholdouts
    for(int nh = 0; nh < nholdouts; ++nh)
    {

//        if strcmp(class(opt.split),'cell')
//            tr = opt.split{nh}.tr;
//            va = opt.split{nh}.va;
//        else
//            tr = opt.split.tr;
//            va = opt.split.va;
//        end
        unsigned long last = lasts[nh];
        copy(tr, indices+nh, indices_cols, 1, indices_mat.rows());
        va = tr+last;
        const unsigned long va_size = indices_cols-last;

//        [n,T]  = size(y(tr,:));

        // here n is last

//        opt.kernel.K = K(tr,tr);
        subK->resize(last, last);
        copy_submatrix(subK->getData(), K.getData(), K.rows(), last, last, tr, tr);


//        opt.predkernel.K = K(va,tr);
        subPredK->resize(va_size, last);
        copy_submatrix(subPredK->getData(), K.getData(), K.rows(), va_size, last, va, tr);

//        opt.predkernel.Ktest = diag(K(va,va));
        T* tmpMat = new T[va_size*va_size];
        subPredKTest->resize(va_size, 1);

        copy_submatrix(tmpMat, K.getData(), K.rows(), va_size, va_size, va, va);
        copy(subPredKTest->getData(), tmpMat, va_size, 1, va_size+1);

        delete[] tmpMat;

        subXtr.resize(last, d);
        subMatrixFromRows(X.getData(), n, d, tr, last, subXtr.getData());

        subYtr.resize(last, t);
        subMatrixFromRows(Y.getData(), n, t, tr, last, subYtr.getData());

        subXva.resize(va_size, d);
        subMatrixFromRows(X.getData(), n, d, va, va_size, subXva.getData());

        subYva.resize(va_size, t);
        subMatrixFromRows(Y.getData(), n, t, va, va_size, subYva.getData());


//        for i = 1:tot
        for(int i=0; i< tot; ++i)
        {
//            opt.paramsel.noises = guesses(i);
            lambda->getData()[0] = guesses[i];

//            opt.rls = rls_gpregr(X(tr,:),y(tr,:),opt);
            GurlsOptionsList* ret_rlsgp = rlsgp.execute(subXtr, subYtr, opt);

            nestedOpt->removeOpt("optimizer");
            nestedOpt->addOpt("optimizer", ret_rlsgp);

//            tmp = pred_gpregr(X(va,:),y(va,:),opt);
            GurlsOptionsList * pred_list = predgp.execute(subXva, subYva, opt);

//            opt.pred = tmp.means;
            nestedOpt->removeOpt("pred");
            nestedOpt->addOpt("pred", pred_list->getOpt("means"));

            pred_list->removeOpt("means", false);

            delete pred_list;


//            opt.perf = opt.hoperf([],y(va,:),opt);
            GurlsOptionsList * perf_list = perfClass->execute(subXva, subYva, opt);
            gMat2D<T>& forho = perf_list->getOptValue<OptMatrix<gMat2D<T> > >("forho");

//            for t = 1:T
//                perf(i,t) = opt.perf.forho(t);
            copy(perf+i, forho.getData(), t, tot, 1);

            delete perf_list;
        }

//        [dummy,idx] = max(perf,[],1);
        work = NULL;
        unsigned long* idx = new unsigned long[t];
        indicesOfMax(perf, tot, t, idx, work, 1);

//        vout.lambdas_round{nh} = guesses(idx);
        T* lambdas_nh = copyLocations(idx, guesses, t, tot);
        copy(lambdas_round + nh, lambdas_nh, t, nholdouts, 1);
        delete [] lambdas_nh;

//        vout.perf{nh} = perf;
        copy(perf_mat->getData()+nh, perf, tot*t, nholdouts, 1);

//        vout.guesses{nh} = guesses;
        copy(ret_guesses+nh, guesses, tot, nholdouts, 1);

    }

    delete nestedOpt;

    delete[] guesses;
    delete perfClass;
    delete[] perf;


    GurlsOptionsList* paramsel;

    if(opt.hasOpt("paramsel"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt("paramsel", opt);

        paramsel = GurlsOptionsList::dynacast(tmp_opt->getOpt("paramsel"));
        tmp_opt->removeOpt("paramsel", false);
        delete tmp_opt;

        paramsel->removeOpt("lambdas_round");
        paramsel->removeOpt("guesses");
        paramsel->removeOpt("perf");
        paramsel->removeOpt("lambdas");
    }
    else
        paramsel = new GurlsOptionsList("paramsel");


    paramsel->addOpt("lambdas_round", new OptMatrix<gMat2D<T> >(*lambdas_round_mat));
    paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*perf_mat));
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));


    gMat2D<T> *l = new gMat2D<T>(1, t);

//    if numel(vout.lambdas_round) > 1
    if(nholdouts>1)
    {
        T *lambdas = new T[t];
//        lambdas = cell2mat(vout.lambdas_round');
//        vout.lambdas = mean(lambdas);
        mean(lambdas_round, lambdas, nholdouts, t, t);

        copy(l->getData(), lambdas, t);

        delete [] lambdas;
    }
    else
    {
//        vout.lambdas = vout.lambdas_round{1};
        copy(l->getData(), lambdas_round, t);
    }

    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*l));

    return paramsel;
}


}

#endif // _GURLS_HOGPREGR_H_
