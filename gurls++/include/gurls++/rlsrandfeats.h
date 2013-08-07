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


#ifndef _GURLS_RLSRANDFEATS_H_
#define _GURLS_RLSRANDFEATS_H_

#include <cmath>

#include "gurls++/optimization.h"

#include "gurls++/gmath.h"
#include "gurls++/gmat2d.h"
#include "gurls++/options.h"
#include "gurls++/optlist.h"
#include "gurls++/utils.h"

namespace gurls {

/**
 * \ingroup Optimization
 * \brief RLSRandFeats is the sub-class of Optimizer that computes a classifier for the primal formulation
 * of RLS using the Random Features approach proposed in:
 *  Ali Rahimi, Ben Recht;
 *  Random Features for Large-Scale Kernel Machines;
 *  in Neural Information Processing Systems (NIPS) 2007.
 * The regularization parameter is set to the one found in opt.paramsel.
 * In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
 */

template <typename T>
class RLSRandFeats: public Optimizer<T>{

public:
    /**
     * Computes a classifier for the primal formulation of RLS using the Random Features approach.
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param struct of options with the following fields:
     *  - paramsel.lambdas (set by the paramsel tasks)
     *  - singlelambda
     *  - randfeats.D
     *  - randfeats.samplesize
     *
     * \return adds to opt the field optimizer, which is a list containing the following fields:
     *  - W: matrix of coefficient vectors of rls estimator for each class
     *  - C: empty matrix
     *  - X: empty matrix
     */
    GurlsOptionsList *execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt);
};


template <typename T>
GurlsOptionsList* RLSRandFeats<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt)
{
    //    lambda = opt.singlelambda(opt.paramsel.lambdas);
    const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
    T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());

    const unsigned long sampleSize = opt.getOptAsNumber("randfeats.samplesize");
    const unsigned long D = opt.getOptAsNumber("randfeats.D");

//    n = size(X,1);
    const unsigned long n = X.rows();
//    const unsigned long d = X.cols();
    const unsigned long t = Y.cols();

//    if or(opt.randfeats.samplesize < 0, opt.randfeats.samplesize > n)
//        ni = n;
//    else
//        ni = opt.randfeats.samplesize;
//    end

    const unsigned long ni = (sampleSize < 0 || sampleSize > n)? n : sampleSize;

    const unsigned long D2 = 2*D;

    T *XtX = new T[D2*D2];
    T *Xty = new T[D2*t];

//    [XtX,Xty,rls.proj] = rp_factorize_large_real(X,y,opt.randfeats.D,ni);
    gMat2D<T> *rls_proj = rp_factorize_large_real(X, Y, D, ni, XtX, Xty);

//    rls.W = rls_primal_driver( XtX, Xty, n, lambda );
    gMat2D<T> *W = rls_primal_driver(XtX, Xty, D2, D2, t, lambda);

    delete [] XtX;
    delete [] Xty;

    GurlsOptionsList *optimizer = new GurlsOptionsList("optimizer");


    optimizer->addOpt("proj", new OptMatrix<gMat2D<T> >(*rls_proj));

//    rls.W = rls_primal_driver( XtX, Xty, n, lambda );
    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

//    rls.C = [];
    gMat2D<T>* emptyC = new gMat2D<T>();
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

//    rls.X = [];
    gMat2D<T>* emptyX = new gMat2D<T>();
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));


    return optimizer;
}

}
#endif // _GURLS_RLSRANDFEATS_H_

