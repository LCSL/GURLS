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


#ifndef _GURLS_PREDKERNELTRAINTEST_H_
#define _GURLS_PREDKERNELTRAINTEST_H_


#include "predkernel.h"
#include "gmath.h"

#include <string>

namespace gurls {


/**
 * \ingroup PredKernels
 * \brief PredKernelTrainTest is the sub-class of PredKernel that computes the kernel matrix between training and test sets
 */

template <typename T>
class PredKernelTrainTest: public PredKernel<T>
{
public:
   /**
     * Computes the kernel matrix between the training set and the test set X, in order to predict the labels for X
     *
     * \param X input data matrix
     * \param Y not used
     * \param opt options with the following required fields:
     *  - optimizer (list with the field X containing the input matrix. Required  only if the field type of kernel is 'rbf' or 'chisquared', is is settable with the class Optimizer and its subclasses )
     *  - kernel (list with the field type, settable through the class Kernel and its subclasses)
     *  - testkernel (only if opt.kernel.type is 'load')
     *  - paramsel (list with the field sigma, required, only if opt.kernel.type is 'rbf', and settable with the class ParamSel and its subclasses SigLam and SiglamHo)
     *
     * \return adds to opt the field predkernel, which is a list with at least the field K containing the kernel matrix
     */

    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
GurlsOptionsList *PredKernelTrainTest<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y, const GurlsOptionsList &opt) throw(gException)
{
    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    const GurlsOptionsList* kernel = GurlsOptionsList::dynacast(opt.getOpt("kernel"));
    const GurlsOptionsList* optimizer = GurlsOptionsList::dynacast(opt.getOpt("optimizer"));

    std::string kernelType = kernel->getOptAsString("type");

    const gMat2D<T>& rls_X_mat = OptMatrix<gMat2D<T> >::dynacast(optimizer->getOpt("X"))->getValue();


    if(X_OMR.cols() != rls_X_mat.cols())
        throw gException(Exception_Inconsistent_Size);


    const int xr = X_OMR.rows();
    const int xc = X_OMR.cols();
    const int rls_xr = rls_X_mat.rows();

    gMat2D<T> rls_X(xc, rls_xr);
    rls_X_mat.transpose(rls_X);


    GurlsOptionsList* predkernel = new GurlsOptionsList("predkernel");
    predkernel->addOpt("type", kernelType);

    gMat2D<T>* K_m;

    if(kernelType == "rbf")
    {
        const GurlsOptionsList* paramsel = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));

        double sigma = paramsel->getOptAsNumber("sigma");

//                opt.predkernel.distance = distance(X',opt.rls.X');
        gMat2D<T> dist(rls_xr, xr);

        T* Xt = new T[xc*xr];
        T* rls_Xt = new T[xc*rls_xr];

        transpose(X.getData(), xr, xc, Xt);
        transpose(rls_X.getData(), rls_xr, xc, rls_Xt);

        distance(Xt, rls_Xt, xc, xr, rls_xr, dist.getData());

        delete[] Xt;
        delete[] rls_Xt;

//                fk.distance = opt.predkernel.distance;
        gMat2D<T>* dist_m = new gMat2D<T>(xr, rls_xr);
        dist.transpose(*dist_m);

        predkernel->addOpt("distance", new OptMatrix<gMat2D<T> > (*dist_m));


//            fk.K = exp(-(opt.predkernel.distance)/(opt.paramsel.sigma^2));
        scal(dist.getSize(), (T)(-1.0/pow(sigma, 2)), dist.getData(), 1);
        exp(dist.getData(), dist.getSize());

        K_m = new gMat2D<T>(xr, rls_xr);
        dist.transpose(*K_m);

        if(!optimizer->hasOpt("L"))
        {
            gMat2D<T> *Ktest = new gMat2D<T>(xr, 1);
            set(Ktest->getData(), (T)1.0, xr);

            predkernel->addOpt("Ktest", new OptMatrix<gMat2D<T> >(*Ktest));
        }

    }

    else if(kernelType == "load")
    {
//            load(opt.testkernel);
        std::string testKernel = opt.getOptAsString("testKernel");

//            fk.K = K_tetr;
        K_m = new gMat2D<T>();
        K_m->load(testKernel);

    }

    else if(kernelType == "chisquared")
    {
        const T epsilon = std::numeric_limits<T>::epsilon();

        gMat2D<T> K(rls_xr, xr);
        T* Kbuf = K.getData();
        set(Kbuf, (T)0.0, K.getSize());

//            for i = 1:size(X,1)
        for(int i=0; i<xr; ++i)
        {
//                for j = 1:size(opt.rls.X,1)
            for(int j=0; j<rls_xr; ++j)
            {

//                    fk.K(i,j) = sum(...
//                                    ( (X(i,:) - opt.rls.X(j,:)).^2 ) ./ ...
//                                    ( 0.5*(X(i,:) + opt.rls.X(j,:)) + eps));
                T sum = 0;
                for(int k=0; k< xc; ++k)
                {
                    const T X_ik = X.getData()[i+(xr*k)];
                    const T rlsX_jk = rls_X.getData()[j+(rls_xr*k)];

                    sum += pow(X_ik - rlsX_jk, 2) / static_cast<T>(((0.5*(X_ik + rlsX_jk)) + epsilon));
                }

                Kbuf[i+(xr*j)] = sum;
            }
        }

        K_m = new gMat2D<T>(xr, rls_xr);
        K.transpose(*K_m);

    }

    else
        throw gException(Exception_Required_Parameter_Missing);

    predkernel->addOpt("K", new OptMatrix<gMat2D<T> > (*K_m));


    return predkernel;
}

}

#endif //_GURLS_PREDKERNELTRAINTEST_H_
