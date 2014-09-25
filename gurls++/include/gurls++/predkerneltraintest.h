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


#ifndef _GURLS_PREDKERNELTRAINTEST_H_
#define _GURLS_PREDKERNELTRAINTEST_H_


#include "gurls++/predkernel.h"
#include "gurls++/gmath.h"
#include "gurls++/utils.h"

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
	///
	/// Default constructor
	///
	PredKernelTrainTest():PredKernel<T>("traintest"){}
	
	///
	/// Clone method
	///
	TaskBase *clone()
	{
		return new PredKernelTrainTest<T>();
	}

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
     * \return predkernel GurlsOptionsList with at least the field K containing the kernel matrix
     */

    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList *PredKernelTrainTest<T>::execute(const gMat2D<T>& X, const gMat2D<T>& /*Y*/, const GurlsOptionsList &opt) throw(gException)
{
    const GurlsOptionsList* optimizer = opt.getOptAs<GurlsOptionsList>("optimizer");

    std::string kernelType = opt.getOptValue<OptString>("kernel.type");

    const gMat2D<T>& rls_X = optimizer->getOptValue<OptMatrix<gMat2D<T> > >("X");


    const unsigned long xr = X.rows();
    const unsigned long xc = X.cols();
    const unsigned long rls_xr = rls_X.rows();

    if(xc != rls_X.cols())
        throw gException(Exception_Inconsistent_Size);


    GurlsOptionsList* predkernel = new GurlsOptionsList("predkernel");
    predkernel->addOpt("type", kernelType);

    gMat2D<T>* K;

    if(kernelType == "rbf")
    {
        double sigma = opt.getOptValue<OptNumber>("paramsel.sigma");

//                opt.predkernel.distance = distance(X',opt.rls.X');
        gMat2D<T> *dist = new gMat2D<T>(xr, rls_xr);

        distance_transposed(X.getData(), rls_X.getData(), xc, xr, rls_xr, dist->getData());


//                fk.distance = opt.predkernel.distance;
        predkernel->addOpt("distance", new OptMatrix<gMat2D<T> > (*dist));


        K = new gMat2D<T>(xr, rls_xr);
        copy(K->getData(), dist->getData(), dist->getSize());

//            fk.K = exp(-(opt.predkernel.distance)/(opt.paramsel.sigma^2));
        scal(K->getSize(), (T)(-1.0/pow(sigma, 2)), K->getData(), 1);
        exp(K->getData(), K->getSize());

        if(optimizer->hasOpt("L"))
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
        K = new gMat2D<T>();
        K->load(testKernel);

    }

    else if(kernelType == "chisquared")
    {
        const T epsilon = std::numeric_limits<T>::epsilon();

        K = new gMat2D<T>(xr, rls_xr);
        T* Kbuf = K->getData();
        set(Kbuf, (T)0.0, K->getSize());

//            for i = 1:size(X,1)
        for(unsigned long i=0; i<xr; ++i)
        {
//                for j = 1:size(opt.rls.X,1)
            for(unsigned long j=0; j<rls_xr; ++j)
            {

//                    fk.K(i,j) = sum(...
//                                    ( (X(i,:) - opt.rls.X(j,:)).^2 ) ./ ...
//                                    ( 0.5*(X(i,:) + opt.rls.X(j,:)) + eps));
                T sum = 0;
                for(unsigned long k=0; k< xc; ++k)
                {
                    const T X_ik = X.getData()[i+(xr*k)];
                    const T rlsX_jk = rls_X.getData()[j+(rls_xr*k)];

                    sum += pow(X_ik - rlsX_jk, 2) / static_cast<T>(((0.5*(X_ik + rlsX_jk)) + epsilon));
                }

                Kbuf[i+(xr*j)] = sum;
            }
        }

    }
    else if(kernelType == "linear")
    {
        //fk.K = X*opt.rls.X';
        K = new gMat2D<T>(xr, rls_xr);
        dot(X.getData(), rls_X.getData(), K->getData(), xr, xc, rls_xr, xc, xr, rls_xr, CblasNoTrans, CblasTrans, CblasColMajor); 

        //if isfield(opt.rls,'L')
        if(optimizer->hasOpt("L"))
        {
            //fk.Ktest = sum(X.^2,2);
            gMat2D<T> *Ktest = new gMat2D<T>(xr, 1);

            sum_col_squared(X.getData(), Ktest->getData(), xr, xc);
            predkernel->addOpt("Ktest", new OptMatrix<gMat2D<T> >(*Ktest));
        }
    }

    else
        throw gException(Exception_Required_Parameter_Missing);

    predkernel->addOpt("K", new OptMatrix<gMat2D<T> > (*K));

    return predkernel;
}

}

#endif //_GURLS_PREDKERNELTRAINTEST_H_
