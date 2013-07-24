/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-1013, IIT@MIT Lab
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


#include "gurls++/norm.h"
#include "gurls++/gmath.h"
#include "gurls++/optmatrix.h"

#include <string>

namespace gurls {

/**
 * \ingroup Norms
 * \brief NormZScore is the sub-class of Norm that centers and rescales the input data matrix X.
 */

template <typename T>
class NormZScore: public Norm<T>
{
public:
    /**
     * Normalizes the input data matrix X, centering them and rescaling it so that each dimension has std. deviation 1. Then saves stats in a file with name root specified in the field name of opt
     * \param X input data matrix
     * \param Y input data matrix
     * \param opt not used
     *
     * \return spheriphied input data matrix with mean and standard deviation
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)  throw(gException);

protected:
    void centerRescale(gMat2D<T> &M, T *stdDevs, const T *means);
};

template<typename T>
GurlsOptionsList* NormZScore<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& /*opt*/) throw(gException)
{
//    [n,d] = size(X);
    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = Y.cols();

//    meanX = mean(X);

    gMat2D<T> *v_meanX = new gMat2D<T>(1, d);
    mean(X.getData(), v_meanX->getData(), n, d, d);

    gMat2D<T> *v_meanY = new gMat2D<T>(1, t);
    mean(Y.getData(), v_meanY->getData(), n, t, t);

//    stdX = std(X) + eps;
//    X = X - repmat(meanX, n, 1);
//    X = X./repmat(stdX, n, 1);

    gMat2D<T> *v_stdX = new gMat2D<T>(1, d);
    gMat2D<T> *v_stdY = new gMat2D<T>(1, t);

    gMat2D<T> *retX = new gMat2D<T>(n, d);
    copy(retX->getData(), X.getData(), retX->getSize());

    gMat2D<T> *retY = new gMat2D<T>(n, t);
    copy(retY->getData(), Y.getData(), retY->getSize());

    centerRescale(*retX, v_stdX->getData(), v_meanX->getData());
    centerRescale(*retY, v_stdY->getData(), v_meanY->getData());

    GurlsOptionsList* norm = new GurlsOptionsList("norm");
    norm->addOpt("X", new OptMatrix<gMat2D<T> >(*retX));
    norm->addOpt("meanX", new OptMatrix<gMat2D<T> >(*v_meanX));
    norm->addOpt("stdX", new OptMatrix<gMat2D<T> >(*v_stdX));

    norm->addOpt("Y", new OptMatrix<gMat2D<T> >(*retY));
    norm->addOpt("meanY", new OptMatrix<gMat2D<T> >(*v_meanY));
    norm->addOpt("stdY", new OptMatrix<gMat2D<T> >(*v_stdY));

    return norm;
}

template<typename T>
void NormZScore<T>::centerRescale(gMat2D<T> &M, T *stdDevs, const T *means)
{
    const unsigned long n = M.rows();
    const unsigned long d = M.cols();

    const T epsilon = std::numeric_limits<T>::epsilon();

    T* column = M.getData();
    T* std_it = stdDevs;
    const T* mean_it = means;
    for(unsigned long i=0; i<d; ++i, column+=n, ++std_it, ++mean_it)
    {
        axpy(n, (T)-1.0, mean_it, 0, column, 1);

        T norm = nrm2(n, column, 1);
        T stdDev = sqrt( (norm*norm) / (n-1)) + epsilon;

        *std_it = stdDev;
        scal(n, (T)1.0/stdDev, column, 1);
    }
}

}

#endif //_GURLS_NORMZSCORE_H_
