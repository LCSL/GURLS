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


#ifndef _GURLS_NORML2_H_
#define _GURLS_NORML2_H_


#include "norm.h"
#include "gmath.h"

#include <limits>

namespace gurls {

/**
 * \ingroup Norms
 * \brief NormL2 is the sub-class of Norm that spheriphies the data according to the l2 norm.
 */

template <typename T>
class NormL2: public Norm<T>
{
public:
    /**
     * Spheriphies the data according to the l2 norm.
     * \param X input data matrix
     * \param Y not used
     * \param opt not used
     * \return spheriphied input data matrix
     */
    gMat2D<T>* execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
gMat2D<T>* NormL2<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& /*Y*/, GurlsOptionsList& /*opt*/) throw(gException)
{
    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    const int m = X_OMR.rows();
    const int n = X_OMR.cols();

    T norm;

//    for j = 1:size(X,1)
    for(int j=0; j<m; ++j)
    {
//        X(j,:) = X(j,:)/(norm(X(j,:)) + eps);
        norm = nrm2(n, X.getData()+j, m) + std::numeric_limits<T>::epsilon();
        scal(n, (T)1.0/norm, X.getData()+j, m);
    }

    gMat2D<T>* ret = new gMat2D<T>(X_OMR.rows(), X_OMR.cols());
    X.transpose(*ret);

    return ret;
}

}

#endif //_GURLS_NORML2_H_
