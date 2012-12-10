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


#ifndef _GURLS_CHISQUAREDKERNEL_H_
#define _GURLS_CHISQUAREDKERNEL_H_


#include "kernel.h"
#include "gmath.h"

namespace gurls {

/**
 * \ingroup Kernels
 * \brief KernelChisquared is the sub-class of Kernel that builds the kernel matrix for a chi-squared model
 */

template <typename T>
class KernelChisquared: public Kernel<T>
{
public:
    /**
     * Builds the symmetric kernel matrix of matrix X for a chi-squared model.
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt not udes
     *
     * \return kernel, a GurslOptionList with the following fields:
     *  - type = "chisquared"
     *  - K = the kernel matrix
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
GurlsOptionsList *KernelChisquared<T>::execute(const gMat2D<T>& X, const gMat2D<T>& /*Y*/, const GurlsOptionsList &/*opt*/) throw(gException)
{
    const int n = X.rows();
    const int t = X.cols();


    gMat2D<T>* K_m = new gMat2D<T>(n, n);
    T* K = K_m->getData();

    const T epsilon = std::numeric_limits<T>::epsilon();

    set(K, (T)0.0, n*n); //diagonal

    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<i; ++j)
        {
            T sum = 0;
            for(int k=0; k< t; ++k)
            {
                const T X_ik = X.getData()[i+(n*k)];
                const T X_jk = X.getData()[j+(n*k)];

                sum += pow(X_ik - X_jk, 2) / static_cast<T>(((0.5*(X_ik + X_jk)) + epsilon));
            }

            K[i+(n*j)] = K[j+(n*i)] = sum;
        }
    }

    //  kernel.type = 'chisquared';
    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");

    kernel->addOpt("type", "chisquared");
    kernel->addOpt("K", new OptMatrix<gMat2D<T> >(*K_m));

    return kernel;
}

}

#endif //_GURLS_CHISQUAREDKERNEL_H_
