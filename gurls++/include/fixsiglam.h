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


#ifndef _GURLS_FIXSIGLAM_H_
#define _GURLS_FIXSIGLAM_H_


#include "options.h"
#include "optlist.h"
#include "gmat2d.h"

#include "paramsel.h"


namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelFixSigLam is the sub-class of ParamSelection that sets the regularization parameters to constants
 */

template <typename T>
class ParamSelFixSigLam: public ParamSelection<T>
{
public:
    /**
     * Sets the regularization parameter lambda  and sigma to the constant 1.
     * \param X not used
     * \param Y not used
     * \param opt not used
     *
     * \return adds the following fields to opt:
     *  - lambdas(=1.0)
     *  - sigma(=1.0)
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList *ParamSelFixSigLam<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");

    gMat2D<T> *lambda = new gMat2D<T>(1,1);
    lambda->getData()[0] = (T)1.0;
    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*lambda));

    paramsel->addOpt("sigma", new OptNumber( 1.0 ));

    return paramsel;
}

}

#endif //_GURLS_FIXSIGLAM_H_
