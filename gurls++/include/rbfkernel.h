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


#ifndef _GURLS_RBFKERNEL_H_
#define _GURLS_RBFKERNEL_H_


#include "kernel.h"
#include "gmath.h"
#include "utils.h"

namespace gurls {

/**
 * \ingroup Kernels
 * \brief KernelRBF is the sub-class of Kernel that builds the Gaussian kernel matrix
 */

template <typename T>
class KernelRBF: public Kernel<T>
{
public:
    /**
     * Builds the symmetric kernel matrix of matrix X for a gaussian model.
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following fields:
     *  - paramsel (list with the required field sigma, settable with the class ParamSelection and its subclasses Siglam and SiglamHo)
     *
     * \return kernel, a GurslOptionList with the following fields:
     *  - type = "rbf"
     *  - K = the kernel matrix
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
GurlsOptionsList *KernelRBF<T>::execute(const gMat2D<T>& X, const gMat2D<T>& /*Y*/, const GurlsOptionsList &opt) throw(gException)
{
    const int xr = X.rows();
    const int xc = X.cols();


//    if ~isfield(opt.kernel,'distance')
//        opt.kernel.distance = distance(X',X');
//        kernel.distance = opt.kernel.distance;
//    end

    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
    gMat2D<T> *dist;

    bool oldDistance = false;

    if(opt.hasOpt("kernel"))
    {
        const GurlsOptionsList* opt_kernel = opt.getOptAs<GurlsOptionsList>("kernel");

        if(opt_kernel->hasOpt("distance"))
            oldDistance = true;
    }

    if(oldDistance)
    {
        const gMat2D<T> &opt_dist = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.distance");

        dist = new gMat2D<T>(opt_dist);
    }
    else
    {
        dist = new gMat2D<T>(xr, xr);

        distance_transposed(X.getData(), X.getData(), xc, xr, xr, dist->getData());
    }


    kernel->addOpt("distance", new OptMatrix<gMat2D<T> >(*dist));

    double sigma = opt.getOptValue<OptNumber>("paramsel.sigma");

    const int len = xr*xr;
    gMat2D<T> *K = new gMat2D<T>(dist->getData(), xr, xr, true);

//    D = -(opt.kernel.distance);
//    K = exp(D/(opt.paramsel.sigma^2));
    scal(len, (T)(-1.0/pow(sigma, 2)),  K->getData(), 1);
    exp(K->getData(), len);

//    kernel.type = 'rbf';
    kernel->addOpt("type", "rbf");

    kernel->addOpt("K", new OptMatrix<gMat2D<T> >(*K));

    return kernel;
}

}

#endif //_GURLS_RBFKERNEL_H_
