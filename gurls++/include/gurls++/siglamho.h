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


#ifndef _GURLS_SIGLAMHO_H_
#define _GURLS_SIGLAMHO_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <set>

#include "gurls++/options.h"
#include "gurls++/optlist.h"
#include "gurls++/gmat2d.h"
#include "gurls++/gvec.h"
#include "gurls++/gmath.h"

#include "gurls++/paramsel.h"
#include "gurls++/perf.h"
#include "gurls++/rbfkernel.h"
#include "gurls++/hodual.h"

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
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList *ParamSelSiglamHo<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //  [n,T]  = size(y);
    const unsigned long t = Y.cols();


    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");
    nestedOpt->copyOpt("nlambda", opt);
    nestedOpt->copyOpt("nholdouts", opt);
    nestedOpt->copyOpt("hoperf", opt);
    nestedOpt->copyOpt("smallnumber", opt);
    nestedOpt->copyOpt("split", opt);

    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
    kernel->addOpt("type", "rbf");
    nestedOpt->addOpt("kernel", kernel);


    GurlsOptionsList* paramsel;

    if(opt.hasOpt("paramsel"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt("paramsel", opt);

        paramsel = GurlsOptionsList::dynacast(tmp_opt->getOpt("paramsel"));
        tmp_opt->removeOpt("paramsel", false);
        delete tmp_opt;

        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("sigma");
    }
    else
        paramsel = new GurlsOptionsList("paramsel");


    gMat2D<T>* dist = NULL;

    // if ~isfield(opt.kernel,'distance')
    if(!kernel->hasOpt("distance"))
        // 	opt.kernel.distance = squareform(pdist(X));
    {
        dist = new gMat2D<T>(X.rows(), X.rows());

        squareform<T>(X.getData(), X.rows(), X.cols(), dist->getData(), X.rows());

        T *distSquared = new T[X.rows()*X.rows()];
        copy(distSquared , dist->getData(), X.rows()*X.rows());

        mult<T>(distSquared, distSquared, dist->getData(), X.rows()*X.rows());

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
        int d_len = X.rows()*(X.rows()-1)/2;
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
        nestedOpt->addOpt("sigmamin", new OptNumber(sqrt( distY[firstPercentile]) ));

        delete [] distY;
    }
    else
    {
        nestedOpt->addOpt("sigmamin", new OptNumber(opt.getOptAsNumber("sigmamin")));
    }

    T sigmamin = static_cast<T>(nestedOpt->getOptAsNumber("sigmamin"));

    //  if ~isfield(opt,'sigmamax')
    if(!opt.hasOpt("sigmamax"))
    {
        // 	%D = sort(opt.kernel.distance);
        // 	%opt.sigmamax = median(D(n,:));
        T mAx = *(std::max_element(dist->getData(),dist->getData()+ dist->getSize()));

        // 	opt.sigmamax = max(max(opt.kernel.distance));
        nestedOpt->addOpt("sigmamax", new OptNumber( sqrt( mAx )));
    }
    else
    {
        nestedOpt->addOpt("sigmamax", new OptNumber(opt.getOptAsNumber("sigmamax")));
    }

    T sigmamax = static_cast<T>(nestedOpt->getOptAsNumber("sigmamax"));

    // if opt.sigmamin <= 0
    if( le(sigmamin, (T)0.0) )
    {
        // 	opt.sigmamin = eps;
        nestedOpt->removeOpt("sigmamin");
        nestedOpt->addOpt("sigmamin", new OptNumber(std::numeric_limits<T>::epsilon()));
        sigmamin = std::numeric_limits<T>::epsilon();
    }

    // if opt.sigmamin <= 0
    if( le(sigmamin, (T)0.0) )
    {
        // 	opt.sigmamax = eps;
        nestedOpt->removeOpt("sigmamax");
        nestedOpt->addOpt("sigmamax", new OptNumber(std::numeric_limits<T>::epsilon()));
        sigmamax = std::numeric_limits<T>::epsilon();
    }


    int nlambda = static_cast<int>(opt.getOptAsNumber("nlambda"));
    int nsigma  = static_cast<int>(opt.getOptAsNumber("nsigma"));

    T q = pow( sigmamax/sigmamin, static_cast<T>(1.0/(nsigma-1.0)));

    // PERF = zeros(opt.nsigma,opt.nlambda,T);
    T* perf = new T[nlambda];

    // sigmas = zeros(1,opt.nsigma);

    KernelRBF<T> rbfkernel;
    ParamSelHoDual<T> hodual;


    T maxPerf = (T)-1.0;
    int m = -1;
    T guess = (T)-1.0;

    T* perf_median = new T[nlambda*t];
    T* guesses_median = new T[nlambda];
    T* row = new T[t];

    const unsigned long nholdouts = static_cast<unsigned long>(opt.getOptAsNumber("nholdouts"));
    T* work = new T[nholdouts];

//    for i = 1:opt.nsigma
    for(int i=0; i<nsigma; ++i)
    {
        nestedOpt->addOpt("paramsel", paramsel);

        paramsel->removeOpt("sigma");
        paramsel->addOpt("sigma", new OptNumber( sigmamin * pow(q, i)));

        // 	opt.kernel = kernel_rbf(X,y,opt);
        GurlsOptionsList* retKernel = rbfkernel.execute(X, Y, *nestedOpt);

        nestedOpt->removeOpt("kernel");
        nestedOpt->addOpt("kernel", retKernel);

        nestedOpt->removeOpt("paramsel", false);

        // 	paramsel = paramsel_hodual(X,y,opt);
        GurlsOptionsList* ret_paramsel = hodual.execute(X, Y, *nestedOpt);


//        PERF(i,:,:) = reshape(median(reshape(cell2mat(paramsel.perf')',opt.nlambda*T,nh),2),T,opt.nlambda)';
        gMat2D<T> &perf_mat = ret_paramsel->getOptValue<OptMatrix<gMat2D<T> > >("perf"); // nholdouts x nlambda*t
        median(perf_mat.getData(), perf_mat.rows(), perf_mat.cols(), 1, perf_median, work);

        for(int j=0;j<nlambda;++j)
        {
            getRow(perf_median, nlambda, t, j, row);
            perf[j] = sumv(row, t);
        }

//        guesses(i,:) = median(cell2mat(paramsel.guesses'),1);
        unsigned long mm = std::max_element(perf, perf + nlambda) - perf;

        if( gt(perf[mm], maxPerf) )
        {
            maxPerf = perf[mm];
            m = i;

            gMat2D<T> &guesses_mat = ret_paramsel->getOptValue<OptMatrix<gMat2D<T> > >("guesses"); // nholdouts x nlambda
            median(guesses_mat.getData(), guesses_mat.rows(), guesses_mat.cols(), 1, guesses_median, work);

            guess = guesses_median[mm];
        }

        delete ret_paramsel;
    }

    delete [] row;
    delete [] work;
    delete [] perf;
    delete [] perf_median;
    delete [] guesses_median;
    delete nestedOpt;


    paramsel->removeOpt("sigma");
    paramsel->addOpt("sigma", new OptNumber( sigmamin * pow(q,m) ));

    // % opt lambda
    // vout.lambdas = guesses(m,n)*ones(1,T);

    gMat2D<T> *LAMBDA = new gMat2D<T>(1, t);
    set(LAMBDA->getData(), guess, t);

    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*LAMBDA));

    return paramsel;

}


}

#endif // _GURLS_SIGLAMHO_H_
