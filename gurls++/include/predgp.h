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

#ifndef _GURLS_PREDGP_H
#define _GURLS_PREDGP_H

#include <cmath>

#include "pred.h"

#include "gmath.h"
#include "gmat2d.h"
#include "options.h"
#include "optlist.h"

using namespace std;

namespace gurls {

/**
 * \ingroup Prediction
 * \brief PredGPRegr is the sub-class of Prediction that computes the predictions of GP
 */

template <typename T>
class PredGPRegr: public Prediction<T> {

public:
    /**
     * computes the predictive distribution of the GP stored in opt.optimizer (L and alpha), on the samples passed in the X matrix.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - Kernel (default)
     *  - L, alpha (settable with the class Optimizers and its subclasses RLSGP)
     *  - predkernel (settable with the class PredKernel and its subclasses PredKernelTrainTest, must contain the fields K and Ktest)
     *
     * \return adds to opt the field pred which contains the following fields:
     *  - means = matrix of output means
     *  - covs = matrix of output covariances
     */
    GurlsOptionsList *execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList *PredGPRegr<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{

//    pred.means = opt.predkernel.K*opt.rls.alpha;

    const GurlsOptionsList* predkernel = GurlsOptionsList::dynacast(opt.getOpt("predkernel"));

    const gMat2D<T> &K_mat = OptMatrix<gMat2D<T> >::dynacast(predkernel->getOpt("K"))->getValue();

    const unsigned long kr = K_mat.rows();
    const unsigned long kc = K_mat.cols();

    gMat2D<T> K(kc, kr);
    K_mat.transpose(K);


    const GurlsOptionsList* rls = GurlsOptionsList::dynacast(opt.getOpt("optimizer"));

    const gMat2D<T> &L_mat = OptMatrix<gMat2D<T> >::dynacast(rls->getOpt("L"))->getValue();
    T* L = new T[L_mat.getSize()];

    const unsigned long lr = L_mat.rows();
    const unsigned long lc = L_mat.cols();

    transpose(L_mat.getData(), lc, lr, L);

    const gMat2D<T> &alpha = OptMatrix<gMat2D<T> >::dynacast(rls->getOpt("alpha"))->getValue();


    gMat2D<T>* means_mat = new gMat2D<T>(kr, alpha.cols());
    dot(K_mat.getData(), alpha.getData(), means_mat->getData(), kr, kc, alpha.rows(), alpha.cols(), kr, alpha.cols(), CblasNoTrans, CblasNoTrans, CblasRowMajor);


    const unsigned long n = Y.rows();

//    pred.vars = zeros(n,1);
    gMat2D<T> *vars_mat = new gMat2D<T>(n, 1);
    T* vars = vars_mat->getData();

    T* v = new T[std::max(kc, n)];

    for(unsigned long i = 0; i<n; ++i)
    {
        getRow(K.getData(), kr, kc, i, v);

////        v = opt.rls.L'\opt.predkernel.K(i,:)';
        mldivide_squared(L, v, lr, lc, kc, 1, CblasTrans);

////        pred.vars(i) = v'*v;
        vars[i] =  dot(kc, v, 1, v, 1);
    }

    delete[] L;

//    pred.vars = opt.predkernel.Ktest - pred.vars;
    const gMat2D<T> &Ktest = OptMatrix<gMat2D<T> >::dynacast(predkernel->getOpt("Ktest"))->getValue();
    copy(v, Ktest.getData(), n);
    axpy(n, (T)-1.0, vars, 1, v, 1);
    copy(vars, v, n);

    delete[] v;

    GurlsOptionsList* pred = new GurlsOptionsList("pred");

    pred->addOpt("means", new OptMatrix<gMat2D<T> >(*means_mat));
    pred->addOpt("vars", new OptMatrix<gMat2D<T> >(*vars_mat));

    return pred;
}

}

#endif // _GURLS_PREDGP_H
