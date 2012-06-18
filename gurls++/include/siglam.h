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


#ifndef _GURLS_SIGLAM_H_
#define _GURLS_SIGLAM_H_

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
#include "loocvdual.h"

namespace gurls {

    /**
     * \brief Siglam is the sub-class of ParamSelection that implements LOO cross-validation with the dual formulation for a rbf kernel
     */

template <typename T>
class Siglam: public ParamSelection<T>{

public:
    /**
     * Performs parameter selection when the dual formulation of RLS is used with rbf kernel.
     * The leave-one-out approach is used oer a 2-dimensional grid of values for the parameters sigma (kernel) and lambda (regularization)
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - nlambda (default)
     *  - nsigma (default)
     *  - hoperf (default)
     *  - smallnumber (default)
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
void Siglam<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
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
        kernel_old = static_cast<GurlsOptionsList*>(opt.getOpt("kernel"));
        opt.removeOpt("kernel", false);
    }

    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
    kernel->addOpt("type", "rbf");
    opt.addOpt("kernel", kernel);

    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");

    gMat2D<T>* dist = new gMat2D<T>(X_OMR.rows(), X_OMR.rows());

    // if ~isfield(opt.kernel,'distance')
    if(!kernel->hasOpt("distance"))
        // 	opt.kernel.distance = squareform(pdist(X));
    {
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

        int firstPercentile = gurls::round( (T)0.01 * d_len + (T)0.5);

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
    if(sigmamin <= 0)
        // 	opt.sigmamin = eps;
        opt.addOpt("sigmamin", new OptNumber(std::numeric_limits<T>::epsilon()));

    // if opt.sigmamin <= 0
    if(sigmamin <= 0)
        // 	opt.sigmamax = eps;
        opt.addOpt("sigmamax", new OptNumber(std::numeric_limits<T>::epsilon()));

    int nlambda = static_cast<int>(opt.getOptAsNumber("nlambda"));
    int nsigma  = static_cast<int>( opt.getOptAsNumber("nsigma"));

    T q = pow( sigmamax/sigmamin, static_cast<T>(1.0/(nsigma-1.0)));

    // LOOSQE = zeros(opt.nsigma,opt.nlambda,T);
    //T* LOOSQE = new T[nsigma*nlambda*t];
    T* acc = new T[nlambda];

    // sigmas = zeros(1,opt.nsigma);
    //     T* sigmas = new T[nsigma];
    // for i = 1:opt.nsigma

    RBFKernel< T > rbfkernel;
    LoocvDual< T > loocvdual;
    T* lwork = new T[nlambda];
    T maxTmp = 0;
    int m = -1;
    T voutLambdas = -1;
    unsigned long* mm = new unsigned long[1];

    for(int i=0;i<nsigma;i++)
    {
        paramsel->removeOpt("sigma",false);
        paramsel->addOpt("sigma", new OptNumber( sigmamin * pow(q, i)));

        // 	opt.kernel = kernel_rbf(X,y,opt);
        opt.addOpt("paramsel",paramsel);
        rbfkernel.execute(X_OMR, Y_OMR, opt);
        //  	GurlsOptionsList* ret_k = static_cast<GurlsOptionsList*>(opt.getOpt("kernel"));
        // 	cout << *ret_k << endl;
        opt.removeOpt("paramsel", false);

        // 	paramsel = paramsel_loocvdual(X,y,opt);
        loocvdual.execute(X_OMR, Y_OMR, opt);


        GurlsOptionsList* ret_paramsel = static_cast<GurlsOptionsList*>(opt.getOpt("paramsel"));


        GurlsOption *looe_opt = ret_paramsel->getOpt("acc");
        GurlsOption *guesses_opt = ret_paramsel->getOpt("guesses");

        gMat2D<T> *looe_mat = &(OptMatrix<gMat2D<T> >::dynacast(looe_opt))->getValue();

        // 	LOOSQE(i,:,:) = paramsel.looe{1};
        // 	guesses(i,:) = paramsel.guesses;
        gMat2D<T> *guesses_mat = &(OptMatrix<gMat2D<T> >::dynacast(guesses_opt))->getValue();

        for(int j=0;j<nlambda;++j)
            acc[j] = sumv<T>(looe_mat->getData() + j*t, t, lwork);

        //        unsigned long* mm = indicesOfMax( acc, t, 1, 1);
        indicesOfMax(acc, nlambda, 1, mm, lwork, 1);

        if( acc[mm[0]] >= maxTmp)
        {
            maxTmp = acc[mm[0]];
            m = i;
            voutLambdas = guesses_mat->operator()(0, mm[0]);
        }

        opt.removeOpt("paramsel");
    }

    delete[] mm;
    delete [] lwork;
    //     delete [] sigmas;
    delete [] acc;

    // M = sum(LOOSQE,3); % sum over classes
    //
    // [dummy,i] = max(M(:));
    // [m,n] = ind2sub(size(M),i);
    //
    // % opt sigma
    // vout.sigma = opt.sigmamin*(q^m);


    paramsel->removeOpt("sigma");
    paramsel->addOpt("sigma", new OptNumber( sigmamin * pow(q,m) ));

    // % opt lambda
    // vout.lambdas = guesses(m,n)*ones(1,T);

    T* voutLam = new T[t];
    set<T>(voutLam,voutLambdas,t);

    OptNumberList* LAMBDA = new OptNumberList();
    for (T* l_it = voutLam, *l_end = voutLam+t; l_it != l_end; ++l_it)
        LAMBDA->add(static_cast<double>(*l_it));

    delete [] voutLam;

    //     opt.addOpt("lambdas", LAMBDA);
    paramsel->addOpt("lambdas", LAMBDA);

    opt.addOpt("paramsel",paramsel);

    opt.removeOpt("kernel");
    if(kernel_old != NULL)
        opt.addOpt("kernel", kernel_old);

}


}

#endif // _GURLS_SIGLAM_H_
