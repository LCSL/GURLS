/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-2013, IIT@MIT Lab
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


#ifndef _GURLS_BIGCALIBRATESGD_H_
#define _GURLS_BIGCALIBRATESGD_H_


#include "bgurls++/bigparamsel.h"
#include "bgurls++/bigarray.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief BigParamselCalibrateSGD is the sub-class of BigParamSelection that implements parameter selection for pegasos
 */

template <typename T>
class BigParamSelCalibrateSGD: public BigParamSelection<T>
{

public:
    /**
     * Performs parameter selection when one wants to solve the problem using rls_pegasos.
     * \param X input data bigarray
     * \param Y labels bigarray
     * \param opt options with the following:
     *  - subsize (default)
     *  - calibfile (default)
     *  - hoperf (default)
     *  - singlelambda (default)
     *  - nlambda (default)
     *
     * \return paramsel, a GurlsOptionList with the following fields:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - W = RLS coefficient vector
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList *BigParamSelCalibrateSGD<T>::execute(const BigArray<T>& /*X*/, const BigArray<T>& /*Y*/, const GurlsOptionsList &/*opt*/)
{
    // TODO

    return new GurlsOptionsList("paramsel");
}


}

#endif // _GURLS_BIGCALIBRATESGD_H_
