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


#ifndef _GURLS_NORMZSCORE_H_
#define _GURLS_NORMZSCORE_H_


#include "norm.h"
#include "gmath.h"

#include <string>

namespace gurls {
    /**
     * \brief NormZScore is the sub-class of Norm that centers and rescales the input data matrix X.
     */
template <typename T>
class NormZScore: public Norm<T>
{
public:
    /**
     * Normalizes the input data matrix X, centering them and rescaling it so that each dimension has std. deviation 1. Then saves stats in a file with name root specified in the field name of opt 
     * \param X input data matrix
     * \param Y not used
     * \param opt not used
     *
     * \return spheriphied input data matrix
     * \return adds the field name to opt
     */
    gMat2D<T>* execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
gMat2D<T>* NormZScore<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& /*Y*/, GurlsOptionsList& opt) throw(gException)
{
    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

//    savevars = {'meanX','stdX'};

//    [n,d] = size(X);
    const int n = X_OMR.rows();
    const int d = X_OMR.cols();

//    meanX = mean(X);

    gMat2D<T> v_meanX (1, d);
    T* meanX = v_meanX.getData();
    mean(X.getData(), meanX, n, d, d);

//    stdX = std(X) + eps;
//    X = X - repmat(meanX, n, 1);
//    X = X./repmat(stdX, n, 1);

    gMat2D<T> v_stdX (1, d);
    T* stdX = v_stdX.getData();

    for(int i=0; i< d; ++i)
    {
        T* column = X.getData()+(n*i);

        axpy(n, (T)-1.0, meanX+i, 0, column, 1);

        stdX[i] = sqrt( pow(nrm2(n, column, 1), 2) / (n-1)) + std::numeric_limits<T>::epsilon();

        scal(n, (T)1.0/stdX[i], column, 1);
    }

//    if numel(savevars) > 0
//        [ST,I] = dbstack();
//        save([opt.name '_' ST(1).name],savevars{:}, '-v7');
//    end

    std::string name = opt.getOptAsString("name");


    std::string fileName;

    fileName = name + "_norm_zscore_meanX.txt";
    v_meanX.save(fileName);

    fileName = name + "_norm_zscore_stdX.txt";
    v_stdX.save(fileName);


    gMat2D<T>* ret = new gMat2D<T>(X_OMR.rows(), X_OMR.cols());
    X.transpose(*ret);

    return ret;

}

}

#endif //_GURLS_NORMZSCORE_H_
