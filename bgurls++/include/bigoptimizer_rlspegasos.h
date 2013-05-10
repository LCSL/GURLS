/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-2013, IIT@MIT Lab
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


#ifndef _GURLS_BIGRLSPEGASOS_H_
#define _GURLS_BIGRLSPEGASOS_H_

#include "bigarray.h"
#include "bigoptimization.h"
#include "utils.h"


namespace gurls {

/**
 * \ingroup Optimization
 * \brief BigRLSPegasos is the sub-class of BigOptimizer that implements the Pegaosos algorithm
 */

template <typename T>
class BigRLSPegasos: public BigOptimizer<T>
{
public:
    /**
     * Computes a classifier for the primal formulation of RLS.
     * The optimization is carried out using a stochastic gradient descent algorithm.
     * The regularization parameter is set to the one found in the field paramsel of opt.
     * In case of multiclass problems, the regularizers need to be combined with the function specified inthe field singlelambda of opt
     *
     * \param X input data bigarray
     * \param Y labels bigarray
     * \param opt options with the following:
     *  - singlelambda (default)
     *  - epochs (default)
     *  - paramsel (settable with the class ParamSelection and its subclasses)
     *  - Xte (test input data matrix, needed for accuracy evaluation)
     *  - yte (test labels matrix, needed for accuracy evaluation)
     *
     * \return returns a list containing the following fields:
     *  - W = matrix of coefficient vectors of rls estimator for each class
     *  - W_sum = sum of the classifiers across iterations
     *  - t0 = stepsize parameter
     *  - count = number of iterations
     *  - acc_last = accuracy of the solution computed in the last iteration
     *  - acc_avg = average accuracy across iterations
     *
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList &opt);
};


template <typename T>
GurlsOptionsList* BigRLSPegasos<T>::execute(const BigArray<T>& /*X*/, const BigArray<T>& /*Y*/, const GurlsOptionsList& /*opt*/)
{
    // TODO

    return new GurlsOptionsList("optimizer");
}

}
#endif // _GURLS_BIGRLSPEGASOS_H_
