/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
  * All rights reserved.
  *
 * author:  M. Santoro
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


#ifndef _GURLS_PARAMSEL_H_
#define _GURLS_PARAMSEL_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "gurls++/options.h"
#include "gurls++/optlist.h"
#include "gurls++/gmat2d.h"
#include "gurls++/gvec.h"
#include "gurls++/gmath.h"

namespace gurls
{

template <typename T>
class ParamSelLoocvPrimal;

template <typename T>
class ParamSelLoocvDual;

template <typename T>
class ParamSelFixLambda;

template <typename T>
class ParamSelFixSigLam;

template <typename T>
class ParamSelCalibrateSGD;

template <typename T>
class ParamSelSiglam;

template <typename T>
class ParamSelSiglamHo;

template <typename T>
class ParamSelHoPrimal;

template <typename T>
class ParamSelHoPrimalr;

template <typename T>
class ParamSelHoDual;

template <typename T>
class ParamSelHoDualr;

template <typename T>
class ParamSelLooGPRegr;

template <typename T>
class ParamSelHoGPRegr;

template <typename T>
class ParamSelSiglamLooGPRegr;

template <typename T>
class ParamSelSiglamHoGPRegr;

/**
 * \ingroup Exceptions
 *
 * \brief BadParamSelectionCreation is thrown when \ref factory tries to generate an unknown parameter selection method
 */
class BadParamSelectionCreation : public std::logic_error
{
public:
    /**
     * Exception constructor.
     */
    BadParamSelectionCreation(std::string type): logic_error("Cannot create type " + type) {}
};

/**
 * \ingroup ParameterSelection
 * \brief ParamSelection is the class that implements parameter selection
 */
template <typename T>
class ParamSelection
{
public:

    /**
     * Implements the selection of the regularization parameter(s)
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the different required fields based on the sub-class
     * \return a GurlsOptionList
     */
    virtual GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) = 0;

    /**
     * Factory function returning a pointer to the newly created object.
     *
     * \warning The returned pointer is a plain, un-managed pointer. The calling
     * function is responsible of deallocating the object.
     */
    static ParamSelection<T>* factory(const std::string& id) throw(BadParamSelectionCreation)
    {
        if(id == "loocvprimal")
            return new ParamSelLoocvPrimal<T>;
        if(id == "loocvdual")
            return new ParamSelLoocvDual<T>;
        if(id == "fixlambda")
            return new ParamSelFixLambda<T>;
        if(id == "calibratesgd")
            return new ParamSelCalibrateSGD<T>;
        if(id == "siglam")
            return new ParamSelSiglam<T>;
        if(id == "siglamho")
            return new ParamSelSiglamHo<T>;
        if(id == "hodual")
            return new ParamSelHoDual<T>;
        if(id == "hodualr")
            return new ParamSelHoDualr<T>;
        if(id == "hoprimal")
            return new ParamSelHoPrimal<T>;
        if(id == "hoprimalr")
            return new ParamSelHoPrimalr<T>;
        if(id == "fixsiglam")
            return new ParamSelFixSigLam<T>;
        if(id == "loogpregr")
            return new ParamSelLooGPRegr<T>;
        if(id == "hogpregr")
            return new ParamSelHoGPRegr<T>;
        if(id == "siglamloogpregr")
            return new ParamSelSiglamLooGPRegr<T>;
        if(id == "siglamhogpregr")
            return new ParamSelSiglamHoGPRegr<T>;

        throw BadParamSelectionCreation(id);
    }
};

}

#endif // _GURLS_PARAMSEL_H_
