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


#ifndef _GURLS_PRIMAL_H
#define _GURLS_PRIMAL_H

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "gmath.h"
#include "options.h"
#include "optlist.h"

#include "pred.h"


namespace gurls {

/**
 * \ingroup Prediction
 * \brief PredPrimal is the sub-class of Prediction that computes the predictions of a linear classifier in the primal formulation
 */

template <typename T>
class PredPrimal: public Prediction<T > {

public:
    /**
     * Computes the predictions of the linear classifier stored in the field optimizer of opt and computed using the primal formulation on the samples passed in the X matrix.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - optimizer (settable with the class Optimizers and its subclasses)
     *
     * \return adds the following fields to opt:
     *  - pred = matrix of predicted labels
     */
   void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void PredPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& /*Y*/,
                                   GurlsOptionsList &opt){
   if (opt.hasOpt("optimizer"))
   {
       GurlsOptionsList* optimizer = GurlsOptionsList::dynacast(opt.getOpt("optimizer"));
//        GurlsOption *g = opt.getOpt("W");
       GurlsOption *g = optimizer->getOpt("W");
       gMat2D<T>& W = OptMatrix< gMat2D<T> >::dynacast(g)->getValue();
       gMat2D<T>* Z = new gMat2D<T>(X.rows(), W.cols());
       *Z = 0;
       dot(X, W, *Z);

       if(opt.hasOpt("pred"))
           opt.removeOpt("pred");

        opt.addOpt("pred", new OptMatrix<gMat2D<T> >(*Z));

   }else {
       throw gException(gurls::Exception_Required_Parameter_Missing);
   }
}


}

#endif // _GURLS_PRIMAL_H
