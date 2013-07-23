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


#ifndef _GURLS_SIGLAMHOGPREGR_H_
#define _GURLS_SIGLAMHOGPREGR_H_

#include <cmath>


#include "gurls++/options.h"
#include "gurls++/optlist.h"
#include "gurls++/gmat2d.h"
#include "gurls++/gvec.h"
#include "gurls++/gmath.h"

#include "gurls++/paramsel.h"
#include "gurls++/hogpregr.h"
#include "gurls++/rbfkernel.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelSiglamHoGPRegr is the sub-class of ParamSelection that implements
 */

template <typename T>
class ParamSelSiglamHoGPRegr: public ParamSelection<T>{

public:
    /**
     *
     */
   GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList* ParamSelSiglamHoGPRegr<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)
{
//    [n,T]  = size(y);
    const unsigned long n = Y.rows();
    const unsigned long t = Y.cols();

    const unsigned long d = X.cols();


    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");
    nestedOpt->copyOpt("nlambda", opt);
    nestedOpt->copyOpt("nholdouts", opt);
    nestedOpt->copyOpt("hoperf", opt);
    nestedOpt->copyOpt("split", opt);
    nestedOpt->copyOpt("singlelambda", opt);


//    if ~isfield(opt,'kernel')
    if(!opt.hasOpt("kernel"))
    {
//        opt.kernel.type = 'rbf';
        GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
        kernel->addOpt("type", "rbf");

        nestedOpt->addOpt("kernel", kernel);
    }
    else
        nestedOpt->copyOpt("kernel", opt);

    GurlsOptionsList* kernel = nestedOpt->getOptAs<GurlsOptionsList>("kernel");

    gMat2D<T> *distance;

//    if ~isfield(opt.kernel,'distance')
    if(!kernel->hasOpt("distance"))
    {
        distance = new gMat2D<T>(n, n);

//        opt.kernel.distance = square_distance(X',X');
        distance_transposed(X.getData(), X.getData(), d, n, n, distance->getData());

        kernel->addOpt("distance", new OptMatrix<gMat2D<T> >(*distance));
    }
    else
        distance = &(kernel->getOptValue<OptMatrix<gMat2D<T> > >("distance"));


//    if ~isfield(opt,'sigmamin')
    if(!opt.hasOpt("sigmamin"))
    {
        int d_len = n*(n-1)/2;
        T* distLinearized = new T[d_len];

//        D = sort(opt.kernel.distance(tril(true(n),-1)));

        const unsigned long size = distance->cols();
        T* it = distLinearized;
        T* d_it = distance->getData();

        for(unsigned long i=1; i< size; ++i)
        {
            gurls::copy(it , d_it+i, size - i);

            it += size - i;
            d_it += size;
        }

        std::sort(distLinearized, distLinearized + d_len);

        //        firstPercentile = round(0.01*numel(D)+0.5);
        int firstPercentile = gurls::round((T)0.01 * d_len + (T)0.5)-1;

        // 	opt.sigmamin = D(firstPercentile);
        nestedOpt->addOpt("sigmamin", new OptNumber(sqrt( distLinearized[firstPercentile]) ));

        delete [] distLinearized;
    }
    else
        nestedOpt->copyOpt("sigmamin", opt);


//    if ~isfield(opt,'sigmamax')
    if(!opt.hasOpt("sigmamax"))
    {
//        opt.sigmamax = sqrt(max(max(opt.kernel.distance)));
        double sigmaMax = sqrt(*std::max_element(distance->getData(), distance->getData()+distance->getSize()));
        nestedOpt->addOpt("sigmamax", new OptNumber(sigmaMax));
    }
    else
        nestedOpt->copyOpt("sigmamax", opt);



//    if opt.sigmamin <= 0
    if( le(nestedOpt->getOptAsNumber("sigmamin"), 0.0))
    {
//        opt.sigmamin = eps;
        nestedOpt->removeOpt("sigmamin");
        nestedOpt->addOpt("sigmamin", new OptNumber(std::numeric_limits<T>::epsilon()));
    }

//    if opt.sigmamax <= 0
    if( le(nestedOpt->getOptAsNumber("sigmamax"), 0.0) )
    {
//        opt.sigmamax = eps;
        nestedOpt->removeOpt("sigmamax");
        nestedOpt->addOpt("sigmamax", new OptNumber(std::numeric_limits<T>::epsilon()));
    }

    const int nsigma = static_cast<int>(opt.getOptAsNumber("nsigma"));
    const int nlambda = static_cast<int>(opt.getOptAsNumber("nlambda"));
    const T sigmamin = static_cast<T>(nestedOpt->getOptAsNumber("sigmamin"));

//    q = (opt.sigmamax/opt.sigmamin)^(1/(opt.nsigma-1));
    T q = static_cast<T>(std::pow(nestedOpt->getOptAsNumber("sigmamax")/nestedOpt->getOptAsNumber("sigmamin"), 1.0/(nsigma-1.0)));

//    PERF = zeros(opt.nsigma,opt.nlambda,T);
    T* perf = new T[nsigma*nlambda];
    set(perf, (T)0.0, nsigma*nlambda);

//    sigmas = zeros(1,opt.nsigma);
//    T* sigmas = new T[nsigma];

    T* guesses = new T[nsigma*nlambda];

    KernelRBF<T> rbf;
    ParamSelHoGPRegr<T> hogp;


    GurlsOptionsList* paramsel_rbf = new GurlsOptionsList("paramsel");
    nestedOpt->addOpt("paramsel", paramsel_rbf);

    T* perf_median = new T[nlambda*t];
    T* guesses_median = new T[nlambda];

    const unsigned long nholdouts = static_cast<unsigned long>(opt.getOptAsNumber("nholdouts"));
    T* work = new T[std::max(nholdouts, t)];

//    for i = 1:opt.nsigma
    for(int i=0; i<nsigma; ++i)
    {
//        sigmas(i) = (opt.sigmamin*(q^(i-1)));
        const T sigma = sigmamin* std::pow(q, i);

//        opt.paramsel.sigma = sigmas(i);
        paramsel_rbf->removeOpt("sigma");
        paramsel_rbf->addOpt("sigma", new OptNumber(sigma));


//        opt.kernel = kernel_rbf(X,y,opt);
        GurlsOptionsList* rbf_kernel = rbf.execute(X, Y, *nestedOpt);

        nestedOpt->removeOpt("kernel");
        nestedOpt->addOpt("kernel", rbf_kernel);


//        paramsel = paramsel_hogpregr(X,y,opt);
        GurlsOptionsList* paramsel_hogp = hogp.execute(X, Y, *nestedOpt);


//        PERF(i,:,:) = reshape(median(reshape(cell2mat(paramsel.perf')',opt.nlambda*T,nh),2),T,opt.nlambda)';
        gMat2D<T> &perf_mat = paramsel_hogp->getOptValue<OptMatrix<gMat2D<T> > >("perf"); //nholdouts x nlambda*t
        median(perf_mat.getData(), perf_mat.rows(), perf_mat.cols(), 1, perf_median, work);

        T* perf_it = perf + i;
        for(int j=0; j<nlambda; ++j, perf_it += nsigma)
        {
            getRow(perf_median, nlambda, t, j, work);

            *perf_it = sumv(work, t);
        }

//        guesses(i,:) = median(cell2mat(paramsel.guesses'),1);
        gMat2D<T> &guesses_mat = paramsel_hogp->getOptValue<OptMatrix<gMat2D<T> > >("guesses"); // nholdouts x nlambda
        median(guesses_mat.getData(), guesses_mat.rows(), guesses_mat.cols(), 1, guesses_median, work);
        copy(guesses+i, guesses_median, nlambda, nsigma, 1);

        delete paramsel_hogp;
    }

    delete [] perf_median;
    delete [] guesses_median;
    delete [] work;

    delete nestedOpt;


//    M = sum(PERF,3); % sum over classes

//    [dummy,i] = max(M(:));
    int i = std::max_element(perf, perf +(nsigma*nlambda)) - perf;

//    [m,n] = ind2sub(size(M),i);
    int im = i%nsigma;

    delete[] perf;

    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");

//    % opt sigma
//    vout.sigma = opt.sigmamin*(q^(m-1));
    paramsel->addOpt("sigma", new OptNumber( sigmamin*(std::pow(q, im))) );

//    % opt lambda
//    vout.noises = guesses(m,n)*ones(1,T);

    gMat2D<T> *lambdas = new gMat2D<T>(1, t);
    set(lambdas->getData(), guesses[i], t);

    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*lambdas));

    delete[] guesses;

    return paramsel;
}


}

#endif // _GURLS_SIGLAMHOGPREGR_H_
