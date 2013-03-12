/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, IIT@MIT Lab
 * All rights reserved.
 *
 * authors:  P.K. Mallapragada, M. Santoro and A. Tacchetti
 * email:   {pavan_m / msantoro / atacchet}@mit.edu
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


#ifndef _GURLS_BIGPREDPRIMAL_H
#define _GURLS_BIGPREDPRIMAL_H

//#include <cstdio>
//#include <cstring>
//#include <iostream>
//#include <cmath>
//#include <algorithm>

//#include "gmath.h"
//#include "options.h"
//#include "optlist.h"

#include "optmatrix.h"

#include "bigarray.h"
#include "bigpred.h"
#include "bigmath.h"

#include <mpi/mpi.h>

namespace gurls {

/**
* \ingroup Prediction
* \brief BigPredPrimal is the sub-class of Prediction that computes the predictions of a linear classifier in the primal formulation
*/

template <typename T>
class BigPredPrimal: public BigPrediction<T >
{

public:
   /**
    * Computes the predictions of the linear classifier stored in the field optimizer of opt and computed using the primal formulation on the samples passed in the X matrix.
    * \param X input data matrix
    * \param Y labels matrix
    * \param opt options with the following:
    *  - optimizer (settable with the class Optimizers and its subclasses)
    *  - tmpfile path of a file used to store and load temporary data
    *  - memlimit maximum amount memory to be used performing matrix multiplications
    *
    * The task can be entirely run in parallel
    *
    * \return a list containing the matrix of predicted labels
    */
  GurlsOptionsList* execute( const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList* BigPredPrimal<T>::execute(const BigArray<T>& X, const BigArray<T>& /*Y*/, const GurlsOptionsList &opt)
{
    const BigArray<T>& W = opt.getOptValue<OptMatrix<BigArray<T> >  >("optimizer.W");

    BigArray<T>* scores = matMult_AB(X, W, opt.getOptAsString("files.pred_filename"), opt.getOptAsNumber("memlimit"));

    GurlsOptionsList* ret = new GurlsOptionsList("pred");
    ret->addOpt("pred", new OptMatrix<BigArray<T> >(*scores));

    return ret;
}


}
#endif // _GURLS_BIGPREDPRIMAL_H
