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


#ifndef _GURLS_RLSGP_H_
#define _GURLS_RLSGP_H_

#include <cmath>

#include "optimization.h"

#include "gmath.h"
#include "gmat2d.h"
#include "options.h"
#include "optlist.h"

namespace gurls {

/**
 * \ingroup Optimization
 * \brief RLSGPRegr is the sub-class of Optimizer that implements GP inference
 */

template <typename T>
class RLSGPRegr: public Optimizer<T>{

public:
    /**
     * Performs GP inference.
     * The noiselevel is set to the one found in the field paramsel of opt.
     * In case of multiclass problems, the noiselevel needs to be combined with the function specified in the field singlelambda of opt
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - singlelambda (default)
     *  - paramsel (settable with the class ParamSelection and its subclasses, and containing field noiselevels)
     *  - kernel (settable with the class Kernel and its subclasses)
     *
     * \return adds to opt the field optimizer, which is a list containing the following fields:
     *  - L
     *  - alpha
     *  - X
     */
    GurlsOptionsList *execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt);
};


template <typename T>
GurlsOptionsList* RLSGPRegr<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)
{
    //    noise = opt.singlelambda(opt.paramsel.lambdas);
    const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
    T noiselevel = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());


    const gMat2D<T> &K_mat = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.K");

    T* K = new T[K_mat.getSize()];
    copy(K, K_mat.getData(), K_mat.getSize());

    //n = size(opt.kernel.K,1);
    const unsigned long n = K_mat.rows();

    //T = size(y,2);
    const unsigned long t = Y.cols();


    //    cfr.L = chol(opt.kernel.K + noise^2*eye(n));
    const T coeff = std::pow(noiselevel, 2);
    unsigned long i=0;
    for(T* it = K; i<n; ++i, it += n+1)
        *it += coeff;

    T* retL = new T[n*n];
    cholesky(K, n, n, retL);

    //    cfr.alpha = cfr.L\(cfr.L'\y);
    gMat2D<T>* alpha = new gMat2D<T>(n, t);
    copy(alpha->getData(), Y.getData(), Y.getSize());

    mldivide_squared(retL, alpha->getData(), n, n, n, t, CblasTrans);
    mldivide_squared(retL, alpha->getData(), n, n, n, t, CblasNoTrans);


    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

//           optimizer.L = L;
    gMat2D<T>* L = new gMat2D<T>(n, n);
    copy(L->getData(), retL, L->getSize());
    optimizer->addOpt("L", new OptMatrix<gMat2D<T> >(*L));

    delete[] retL;

//           optimizer.alpha = alpha;
    optimizer->addOpt("alpha", new OptMatrix<gMat2D<T> >(*alpha));

//    cfr.X = X;
    gMat2D<T>* optX = new gMat2D<T>(X);
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*optX));

    return optimizer;
}

}
#endif // _GURLS_RLSGP_H_

