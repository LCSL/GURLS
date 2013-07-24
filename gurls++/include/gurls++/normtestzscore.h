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


#ifndef _GURLS_NORMTESTZSCORE_H_
#define _GURLS_NORMTESTZSCORE_H_


#include "gurls++/norm.h"
#include "gurls++/gmath.h"
#include "gurls++/optmatrix.h"

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
     * Spheriphies test data using the same mean and std deviation computed for the training set,
     * previously computed and stored in opt by the NormZScore sub-class of the class Norm.
     * \param X input data matrix
     * \param Y input data matrix
     * \param opt option list with the following field:
     *  - meanX
     *  - meanY
     *  - stdX
     *  - stdY
     *
     * \return spheriphied input data matrix
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)  throw(gException);

protected:
    void centerRescale(gMat2D<T> &M, const T *stdDevs, const T *means);
};

template<typename T>
GurlsOptionsList* NormTestZScore<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) throw(gException)
{
//    [n,d] = size(X);
    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long m = Y.rows();
    const unsigned long t = Y.cols();

    GurlsOptionsList* norm = new GurlsOptionsList("norm");

    if(n > 0ul && d > 0ul)
    {
        const gMat2D<T> &v_meanX = opt.getOptValue<OptMatrix<gMat2D<T> > >("meanX");
        const gMat2D<T> &v_stdX = opt.getOptValue<OptMatrix<gMat2D<T> > >("stdX");

        gMat2D<T>* retX = new gMat2D<T>(n, d);
        copy(retX->getData(), X.getData(), retX->getSize());

        centerRescale(*retX, v_stdX.getData(), v_meanX.getData());
        norm->addOpt("X", new OptMatrix<gMat2D<T> >(*retX));
    }

    if(m > 0ul && t > 0ul)
    {
        const gMat2D<T> &v_meanY = opt.getOptValue<OptMatrix<gMat2D<T> > >("meanY");
        const gMat2D<T> &v_stdY = opt.getOptValue<OptMatrix<gMat2D<T> > >("stdY");

        gMat2D<T>* retY = new gMat2D<T>(m, t);
        copy(retY->getData(), Y.getData(), retY->getSize());

        centerRescale(*retY, v_stdY.getData(), v_meanY.getData());
        norm->addOpt("Y", new OptMatrix<gMat2D<T> >(*retY));
    }

    return norm;
}

template<typename T>
void NormTestZScore<T>::centerRescale(gMat2D<T> &M, const T *stdDevs, const T *means)
{
    const unsigned long n = M.rows();
    const unsigned long d = M.cols();

    //    X = X - repmat(meanX, n, 1);
    //    X = X./repmat(stdX, n, 1);
    T* column = M.getData();
    const T* std_it = stdDevs;
    const T* mean_it = means;
    for(unsigned long i=0; i<d; ++i, column+=n, ++std_it, ++mean_it)
    {
        axpy(n, (T)-1.0, mean_it, 0, column, 1);
        scal(n, (T)1.0/(*std_it), column, 1);
    }
}

}

#endif //_GURLS_NORMTESTZSCORE_H_
