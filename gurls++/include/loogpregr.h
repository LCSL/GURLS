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


#ifndef _GURLS_LOOGPREGR_H_
#define _GURLS_LOOGPREGR_H_

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
 * \brief ParamSelLooGPRegr is the sub-class of ParamSelection that implements
 */

template <typename T>
class ParamSelLooGPRegr: public ParamSelection<T>{

public:
    /**
     *
     */
   void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void ParamSelLooGPRegr<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{

    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);

//    [n,T]  = size(y);
    const unsigned long n = Y_OMR.rows();
    const unsigned long t = Y_OMR.cols();

    const unsigned long d = X_OMR.cols();

//    tot = opt.nlambda;
    int tot = static_cast<int>(opt.getOptAsNumber("nlambda"));

//    K = opt.kernel.K;
    GurlsOptionsList *kernel = GurlsOptionsList::dynacast(opt.getOpt("kernel"));
    gMat2D<T> &K_mat = OptMatrix<gMat2D<T> >::dynacast(kernel->getOpt("K"))->getValue();

    gMat2D<T> K(K_mat.cols(), K_mat.rows());
    K_mat.transpose(K);


//    lmax = mean(std(y));

    T* work = new T[t+n];
    T* stdY = new T[t];

    stdDev(Y.getData(), n, t, stdY, work);

    const T lmax = sumv(stdY, t, work)/((T)t);

    delete[] work;
    delete[] stdY;

//    lmin = mean(std(y))*10^-5;
    const T lmin = lmax * (T)1.0e-5;

//    guesses = lmin.*(lmax/lmin).^linspace(0,1,tot);
    T* guesses = new T[tot];
    T* linspc = new T[tot];
    linspace((T)0.0, (T)1.0, tot, linspc);
    const T coeff = lmax/lmin;

    for(int i=0; i< tot; ++i)
        guesses[i] = lmin* std::pow(coeff, linspc[i]);

    delete[] linspc;


//    perf = zeros(tot,T);
    T* perf = new T[tot*t];
    set(perf, (T)0.0, tot*t);

    const int tr_size = n-1;

    unsigned long* tr = new unsigned long[tr_size+1]; // + 1 cell for convenience
    unsigned long* tr_it = tr;
    for(unsigned long i=1; i< n; ++i, ++tr_it)
        *tr_it = i;

    opt.removeOpt("kernel", false);
    GurlsOption* predKernel = NULL;
    if(opt.hasOpt("predkernel"))
    {
        predKernel = opt.getOpt("predkernel");
        opt.removeOpt("predkernel", false);
    }

    GurlsOption* optimizer = NULL;
    if(opt.hasOpt("optimizer"))
    {
        optimizer = opt.getOpt("optimizer");
        opt.removeOpt("optimizer", false);
    }

    gMat2D<T>* tmpK = new gMat2D<T>(tr_size, tr_size);
    gMat2D<T>* tmpPredK = new gMat2D<T>(1, tr_size);
    gMat2D<T>* tmpPredKTest = new gMat2D<T>(1, 1);

    GurlsOptionsList* tmpPredKernel = new GurlsOptionsList("predkernel");
    GurlsOptionsList* tmpKernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* tmpParamSel = new GurlsOptionsList("paramsel");

    opt.addOpt("kernel", tmpKernel);
    opt.addOpt("predkernel", tmpPredKernel);
    opt.addOpt("paramsel", tmpParamSel);

    tmpKernel->addOpt("K", new OptMatrix<gMat2D<T> > (*tmpK));
    tmpPredKernel->addOpt("K", new OptMatrix<gMat2D<T> > (*tmpPredK));
    tmpPredKernel->addOpt("Ktest", new OptMatrix<gMat2D<T> > (*tmpPredKTest));

    gMat2D<T> rlsX(tr_size, d);
    gMat2D<T> rlsY(tr_size, t);

    T* tmpMat = new T[ tr_size * std::max(d, t)];

    gMat2D<T> predX(1, d);
    gMat2D<T> predY(1, t);

    RLSGPRegr<T> rlsgp;
    PredGPRegr<T> predgp;
    Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));


//    for k = 1:n;
    for(unsigned long k = 0; k<n; ++k)
    {
//        tr = setdiff(1:n,k);

//        opt.kernel.K = K(tr,tr);
        copy_submatrix(tmpK->getData(), K.getData(), K_mat.rows(), tr_size, tr_size, tr, tr);

//        opt.predkernel.K = K(k,tr);
        copy_submatrix(tmpPredK->getData(), K.getData(), K_mat.rows(), 1, tr_size, &k , tr);

//        opt.predkernel.Ktest = K(k,k);
        tmpPredKTest->getData()[0] = K.getData()[(k*K_mat.rows()) + k];

//        for i = 1:tot
        for(int i=0; i< tot; ++i)
        {
//            opt.paramsel.noises = guesses(i);
            tmpParamSel->removeOpt("lambdas");
            tmpParamSel->addOpt("lambdas", new OptNumberList(guesses[i]));

//            opt.rls = rls_gpregr(X(tr,:),y(tr,:),opt);
            subMatrixFromRows(X.getData(), n, d, tr, tr_size, tmpMat);
            transpose(tmpMat, tr_size, d, rlsX.getData());

            subMatrixFromRows(Y.getData(), n, t, tr, tr_size, tmpMat);
            transpose(tmpMat, tr_size, t, rlsY.getData());

            rlsgp.execute(rlsX, rlsY, opt);

//            tmp = pred_gpregr(X(k,:),y(k,:),opt);
            getRow(X.getData(), n, d, k, predX.getData());
            getRow(Y.getData(), n, t, k, predY.getData());

            predgp.execute(predX, predY, opt);

//            opt.pred = tmp.means;
            opt.removeOpt("optimizer");
            GurlsOptionsList * pred_list = GurlsOptionsList::dynacast(opt.getOpt("pred"));
            opt.removeOpt("pred", false);

            opt.addOpt("pred", pred_list->getOpt("means"));
            pred_list->removeOpt("means", false);
            delete pred_list;

//            opt.perf = opt.hoperf([],y(k,:),opt);
            perfClass->execute(predX, predY, opt);

            opt.removeOpt("pred");

            GurlsOptionsList * perf_list = GurlsOptionsList::dynacast(opt.getOpt("perf"));
            gMat2D<T>& forho = OptMatrix<gMat2D<T> >::dynacast(perf_list->getOpt("forho"))->getValue();

//            for t = 1:T
            for(unsigned long j = 0; j<t; ++j)
//                perf(i,t) = opt.perf.forho(t)+perf(i,t)./n;
                perf[i+(tot*j)] += forho.getData()[j]/n;


            opt.removeOpt("perf");
        }

        tr[k] = k;

    }

    delete perfClass;

    opt.removeOpt("kernel");
    opt.removeOpt("predkernel");
    opt.removeOpt("paramsel");
//    opt.removeOpt("optimizer");

    opt.addOpt("kernel", kernel);
    if(predKernel != NULL)
        opt.addOpt("predkernel", predKernel);

    if(optimizer != NULL)
        opt.addOpt("optimizer", optimizer);


    GurlsOptionsList* paramsel = NULL;
    if(!opt.hasOpt("paramsel"))
    {
        paramsel = new GurlsOptionsList("paramsel");
        opt.addOpt("paramsel", paramsel);
    }
    else
    {
        paramsel = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));

        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("perf");
        paramsel->removeOpt("guesses");

    }

//    [dummy,idx] = max(perf,[],1);
    unsigned long* idx = new unsigned long[t];
    work = NULL;
    indicesOfMax(perf, tot, t, idx, work, 1);


//    vout.noises = 	guesses(idx);
    T* noises = copyLocations(idx, guesses, t, tot);

    delete[] idx;

    OptNumberList* lambdas = new OptNumberList();
    for (T* l_it = noises, *l_end = noises+t; l_it != l_end; ++l_it)
        lambdas->add(static_cast<double>(*l_it));

    paramsel->addOpt("lambdas", lambdas);

    delete[] noises;

//    vout.perf = 	perf;
    gMat2D<T> *perf_mat = new gMat2D<T>(tot, t);
    transpose(perf, tot, t, perf_mat->getData());
    delete[] perf;

    paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*perf_mat));


//    vout.guesses = guesses;
    gMat2D<T> *guesses_mat = new gMat2D<T>(guesses, tot, 1, true);
    delete[] guesses;

    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

}


}

#endif // _GURLS_LOOGPREGR_H_
