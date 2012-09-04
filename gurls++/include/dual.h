/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, Matteo Santoro
 * All rights reserved.
 *
 * author:  M. Santoro
 * email:   matteo.santoro@gmail.com
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


#ifndef _GURLS_DUAL_H
#define _GURLS_DUAL_H

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "gmath.h"
#include "options.h"
#include "optlist.h"

#include "pred.h"
#include "primal.h"


namespace gurls {

/**
 * \ingroup Prediction
 * \brief PredDual is the sub-class of Prediction that computes the predictions of a linear classifier in the dual formulation
 */

template <typename T>
class PredDual: public Prediction<T> {

public:
    /**
     * Computes the predictions of the linear classifier stored in opt.W and computed using the primal formulation, on the samples passed in the X matrix.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - Kernel (default)
     *  - C, X, W (settable with the class Optimizers and its subclasses RLSDual)
     *  - predkernel (required only if the subfield type of Kernel is different than "linear", and settable with the class PredKernel and its subclasses PredKernelTrainTest)
     *
     * \return adds the following field to opt:
     *  - pred = matrix of predicted labels
     */
    OptMatrix<gMat2D<T> >* execute( const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
OptMatrix<gMat2D<T> >* PredDual<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)
{
    if(opt.hasOpt("kernel"))
    {
        if(opt.getOptValue<OptString>("kernel.type") == "linear")
        {
            PredPrimal<T> pred;
            return pred.execute(X, Y, opt);
        }
    }

    const gMat2D<T> &K = opt.getOptValue<OptMatrix<gMat2D<T> > >("predkernel.K");
    const gMat2D<T> &C = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.C");

    gMat2D<T>* Z = new gMat2D<T>(K.rows(), C.cols());


    dot(K.getData(), C.getData(), Z->getData(), K.rows(), K.cols(), C.rows(), C.cols(), Z->rows(), Z->cols(), CblasNoTrans, CblasNoTrans, CblasRowMajor);

    return new OptMatrix<gMat2D<T> >(*Z);
}

}

#endif // _GURLS_DUAL_H
