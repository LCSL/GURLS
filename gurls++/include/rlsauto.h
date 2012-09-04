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


#ifndef _GURLS_RLSAUTO_H_
#define _GURLS_RLSAUTO_H_

#include "optimization.h"


namespace gurls {

/**
 * \ingroup Optimization
 * \brief RLSAuto is the sub-class of Optimizer that implements automatic selection of primal/dual procedure.
 */

template <typename T>
class RLSAuto: public Optimizer<T>{

public:
    /**
      * computes a RLS classifier, with automatic selection of primal/dual procedure.
      * The regularization parameter is set to the one found in the field paramsel of opt
     * In case of multiclass problems, the regularizers need to be combined with the function specified inthe field singlelambda of opt
      *
      * \param X input data matrix
      * \param Y labels matrix
      * \param opt options with the following:
      *  - singlelambda (default)
      *  - paramsel (settable with the class ParamSelection and its subclasses)
      *  - kernel (list with the filed type set to 'linear', settable with the class Kernel and its KernelLinear)
      *
     * \return adds to opt the field optimizer, which is a list containing the following fields:
      *  - W = matrix of coefficient vectors of rls estimator for each class
      *  - C = empty matrix
      *  - X = empty matrix
      */
    GurlsOptionsList *execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt);
};

template <typename T>
GurlsOptionsList* RLSAuto<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, const GurlsOptionsList& opt)
{
    //[n,d] = size(X);
    const unsigned long n = X_OMR.rows();
    const unsigned long d = X_OMR.cols();

//     if (n > d) % Do primal
    if(n > d)
    {
        //    cfr = rls_primal(X, y, opt);
        RLSPrimal<T> rlsprimal;
        return rlsprimal.execute(X_OMR, Y_OMR, opt);
    }
//  else % Do dual
    else
    {
        //      cfr = rls_dual(X, y, opt);
        RLSDual<T> rlsdual;
        return rlsdual.execute(X_OMR, Y_OMR, opt);
    }
}

}
#endif // _GURLS_RLSAUTO_H_

