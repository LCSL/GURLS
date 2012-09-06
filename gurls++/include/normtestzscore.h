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


#ifndef _GURLS_NORMTESTZSCORE_H_
#define _GURLS_NORMTESTZSCORE_H_


#include "norm.h"
#include "gmath.h"

namespace gurls {

/**
 * \ingroup Norms
 * \brief NormTestZScore is the sub-class of Norm that spheriphies the data according to the statistics cmoputed on the training set.
 */

template <typename T>
class NormTestZScore: public Norm<T>
{
public:
    /**
     * Spheriphies test data using the same mean and std deviation computed for the training set, previously computed and stored in the file with name root specified in the field name of opt by the NormZScore sub-class of the class Norm.
     * \param X input data matrix
     * \param Y not used
     * \param opt option list with the following field:
     *  - name (settable with the NormZScore sub-class of the class Norm)
     *
     * \return spheriphied input data matrix
     */
    gMat2D<T>* execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
gMat2D<T>* NormTestZScore<T>::execute(const gMat2D<T>& X, const gMat2D<T>& /*Y*/, GurlsOptionsList& opt) throw(gException)
{
//    [n,d] = size(X);
    const unsigned long n = X.rows();
    const unsigned long d = X.cols();

    std::string name = opt.getOptAsString("name");

    gMat2D<T> v_meanX;
    gMat2D<T> v_stdX;

    std::string fileName;

    fileName = name + "_norm_zscore_meanX.txt";
    v_meanX.load(fileName);

    fileName = name + "_norm_zscore_stdX.txt";
    v_stdX.load(fileName);

    T* meanX = v_meanX.getData();
    T* stdX = v_stdX.getData();

    gMat2D<T>* retX = new gMat2D<T>(n, d);
    copy(retX->getData(), X.getData(), retX->getSize());

//    X = X - repmat(meanX, n, 1);
//    X = X./repmat(stdX, n, 1);
    for(unsigned long i=0; i<d; ++i)
    {
        T* column = retX->getData()+(n*i);

        axpy(n, (T)-1.0, meanX+i, 0, column, 1);
        scal(n, (T)1.0/stdX[i], column, 1);
    }

    return retX;
}

}

#endif //_GURLS_NORMTESTZSCORE_H_
