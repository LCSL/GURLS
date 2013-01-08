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


#ifndef _GURLS_BIGHOPRIMAL_H_
#define _GURLS_BIGHOPRIMAL_H_


#include "bigparamsel.h"


namespace gurls
{

/**
 * \ingroup ParameterSelection
 * \brief BigParamSelHoPrimal is the subclass of BigParamSelection that implements hold-out cross validation with the primal formulation of RLS
 */
template <typename T>
class BigParamSelHoPrimal: public BigParamSelection<T>
{

public:

    /**
     * Performs parameter selection when the primal formulation of RLS is used.
     * The hold-out approach is used.
     * The performance measure specified by opt.hoperf is maximized.
     * \param X input data bigarray
     * \param Y labels bigarray
     * \param opt options with the following:
     *  - nlambda (default)
     *  - hoperf (default)
     *  - smallnumber (default)
     *  - split (settable with the class Split and its subclasses)
     *
     * \return a GurlsOptionList with the following fields:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - guesses = array of guesses for the regularization parameter lambda
     *  - forho = matrix of validation accuracies for each lambda guess and for each class
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt);

};

template <typename T>
GurlsOptionsList *BigParamSelHoPrimal<T>::execute(const BigArray<T> &/*X*/, const BigArray<T> &/*Y*/, const GurlsOptionsList &/*opt*/)
{
    // TODO

    return new GurlsOptionsList("paramsel");
}

}

#endif // _GURLS_BIGHOPRIMAL_H_
