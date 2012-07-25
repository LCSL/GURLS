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


#ifndef _GURLS_SIGLAMHO_H_
#define _GURLS_SIGLAMHO_H_

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
#include "rbfkernel.h"
#include "hodual.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelSiglam is the sub-class of ParamSelection that implements  hold-out cross validation with the dual formulation for a rbf kernel
 */

template <typename T>
class ParamSelSiglamHo: public ParamSelection<T>{

public:

    /**
     * Performs parameter selection when the dual formulation of RLS is used with rbf kernel.
     * The hold-out approach is used oer a 2-dimensional grid of values for the parameters sigma (kernel) and lambda (regularization)
     * The performance measure specified by opt.hoperf is maximized.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - nlambda (default)
     *  - nsigma (default)
     *  - hoperf (default)
     *  - smallnumber (default)
     *  - split (settable with the class Split and its subclasses)
     *
     * \return adds the field paramsel to opt, which is alist containing the following fields:
     *  - lambdas = array containing the value of the regularization parameter lambda maximizing the mean validation accuracy over the classes, replicated as many times as the number of classes
     *  - sigma = values of the kernel parameter maximizing the validation accuracy
     *  - guesses = array of guesses for the regularization parameter lambda
     *  - acc = matrix of validation accuracies for each lambda guess and for each class
     */
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void ParamSelSiglamHo<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{
    //  [n,T]  = size(y);
//    const int n = Y_OMR.rows();
    const int t = Y_OMR.cols();

    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);

    GurlsOptionsList* kernel_old = NULL;
    if (opt.hasOpt("kernel"))
    {
        kernel_old = GurlsOptionsList::dynacast(opt.getOpt("kernel"));
        opt.removeOpt("kernel", false);
    }

    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
    kernel->addOpt("type", "rbf");
    opt.addOpt("kernel", kernel);


//    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");
    GurlsOptionsList* paramsel = NULL;
    if(!opt.hasOpt("paramsel"))
        paramsel = new GurlsOptionsList("paramsel");
    else
    {
        paramsel = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));

        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("sigma");
    }

    gMat2D<T>* dist = NULL;

    // if ~isfield(opt.kernel,'distance')
    if(!kernel->hasOpt("distance"))
        // 	opt.kernel.distance = squareform(pdist(X));
    {
        dist = new gMat2D<T>(X_OMR.rows(), X_OMR.rows());

        squareform<T>(X.getData(), X_OMR.rows(), X_OMR.cols(), dist->getData(), X_OMR.rows());

        T *distSquared = new T[X_OMR.rows()*X_OMR.rows()];
        copy(distSquared , dist->getData(), X_OMR.rows()*X_OMR.rows());

        mult<T>(distSquared, distSquared, dist->getData(), X_OMR.rows()*X_OMR.rows());

        kernel->addOpt("distance", new OptMatrix<gMat2D<T> >(*dist));
        delete [] distSquared;

    }
    else
    {
        GurlsOption *dist_opt = kernel->getOpt("distance");
        dist = &(OptMatrix<gMat2D<T> >::dynacast(dist_opt))->getValue();
    }


    //  if ~isfield(opt,'sigmamin')
    if(!opt.hasOpt("sigmamin"))
    {
        // 	%D = sort(opt.kernel.distance);
        // 	%opt.sigmamin = median(D(2,:));
        // 	D = sort(squareform(opt.kernel.distance));
        int d_len = X_OMR.rows()*(X_OMR.rows()-1)/2;
        T* distY = new T[d_len];
        //        squareform<T>(dist->getData(), dist->rows(), dist->cols(), distY, 1);

        const int size = dist->cols();
        T* it = distY;

        for(int i=1; i< size; it+=i, ++i)
            copy(it , dist->getData()+(i*size), i);

        std::sort(distY, distY + d_len);
        // 	firstPercentile = round(0.01*numel(D)+0.5);

        int firstPercentile = gurls::round((T)0.01 * d_len + (T)0.5)-1;

        // 	opt.sigmamin = D(firstPercentile);
        opt.addOpt("sigmamin", new OptNumber(sqrt( distY[firstPercentile]) ));

        delete [] distY;
    }

    T sigmamin = static_cast<T>(opt.getOptAsNumber("sigmamin"));

    //  if ~isfield(opt,'sigmamax')
    if(!opt.hasOpt("sigmamax"))
    {
        // 	%D = sort(opt.kernel.distance);
        // 	%opt.sigmamax = median(D(n,:));
        T mAx = *(std::max_element(dist->getData(),dist->getData()+ dist->getSize()));

        // 	opt.sigmamax = max(max(opt.kernel.distance));
        opt.addOpt("sigmamax", new OptNumber( sqrt( mAx )));
    }

    T sigmamax = static_cast<T>(opt.getOptAsNumber("sigmamax"));

    // if opt.sigmamin <= 0
    if( le(sigmamin, (T)0.0) )
        // 	opt.sigmamin = eps;
        opt.addOpt("sigmamin", new OptNumber(std::numeric_limits<T>::epsilon()));

    // if opt.sigmamin <= 0
    if( le(sigmamin, (T)0.0) )
        // 	opt.sigmamax = eps;
        opt.addOpt("sigmamax", new OptNumber(std::numeric_limits<T>::epsilon()));


    int nlambda = static_cast<int>(opt.getOptAsNumber("nlambda"));
    int nsigma  = static_cast<int>(opt.getOptAsNumber("nsigma"));

    T q = pow( sigmamax/sigmamin, static_cast<T>(1.0/(nsigma-1.0)));

    // PERF = zeros(opt.nsigma,opt.nlambda,T);
    T* perf = new T[nlambda];

    // sigmas = zeros(1,opt.nsigma);
//     T* sigmas = new T[nsigma];

    KernelRBF<T> rbfkernel;
    ParamSelHoDual<T> hodual;


    T maxPerf = (T)-1.0;
    int m = -1;
    T guess = (T)-1.0;

    T* perf_median = new T[nlambda*t];
    T* guesses_median = new T[nlambda];
    T* row = new T[t];

    const int nholdouts = static_cast<int>(opt.getOptAsNumber("nholdouts"));
    T* work = new T[std::max(nholdouts, t+1)];

//    for i = 1:opt.nsigma
    for(int i=0; i<nsigma; ++i)
    {
        paramsel->removeOpt("sigma");
        paramsel->addOpt("sigma", new OptNumber( sigmamin * pow(q, i)));

        // 	opt.kernel = kernel_rbf(X,y,opt);
        opt.addOpt("paramsel", paramsel);
        rbfkernel.execute(X_OMR, Y_OMR, opt);

        opt.removeOpt("paramsel", false);

        // 	paramsel = paramsel_loocvdual(X,y,opt);
        hodual.execute(X_OMR, Y_OMR, opt);


        GurlsOptionsList* ret_paramsel = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));

//        PERF(i,:,:) = reshape(median(reshape(cell2mat(paramsel.perf')',opt.nlambda*T,nh),2),T,opt.nlambda)';
        gMat2D<T> &perf_mat = (OptMatrix<gMat2D<T> >::dynacast(ret_paramsel->getOpt("perf")))->getValue(); // nholdouts x nlambda*t
        median(perf_mat.getData(), perf_mat.cols(), perf_mat.rows(), 2, perf_median, work); //inverted parameters because perf_mat is in row-major order

        for(int j=0;j<nlambda;++j)
        {
            getRow(perf_median, nlambda, t, j, row);
            perf[j] = sumv(row, t, work);
        }

//        guesses(i,:) = median(cell2mat(paramsel.guesses'),1);
        unsigned long mm = std::max_element(perf, perf + nlambda) - perf;

        if( gt(perf[mm], maxPerf) )
        {
            maxPerf = perf[mm];
            m = i;

            gMat2D<T> &guesses_mat = (OptMatrix<gMat2D<T> >::dynacast(ret_paramsel->getOpt("guesses")))->getValue(); // nholdouts x nlambda
            median(guesses_mat.getData(), guesses_mat.cols(), guesses_mat.rows(), 2, guesses_median, work); //inverted parameters because guesses_mat is in row-major order

            guess = guesses_median[mm];
        }

        opt.removeOpt("paramsel");
    }

    delete [] row;
    delete[] work;
//    delete [] sigmas;
    delete [] perf;
    delete [] perf_median;
    delete [] guesses_median;


    paramsel->removeOpt("sigma");
    paramsel->addOpt("sigma", new OptNumber( sigmamin * pow(q,m) ));

    // % opt lambda
    // vout.lambdas = guesses(m,n)*ones(1,T);

    OptNumberList* LAMBDA = new OptNumberList();
    const double lambda = static_cast<double>(guess);
    for (int i=0; i<t; ++i)
        LAMBDA->add(lambda);

    paramsel->addOpt("lambdas", LAMBDA);

    opt.addOpt("paramsel",paramsel);

    opt.removeOpt("kernel");
    if(kernel_old != NULL)
        opt.addOpt("kernel", kernel_old);

}


}

#endif // _GURLS_SIGLAMHO_H_
