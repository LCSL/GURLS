/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011, IIT@MIT Lab
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

#include "options.h"
#include "optlist.h"
#include "gmat2d.h"
#include "gvec.h"
#include "gmath.h"

namespace gurls {

template <typename T>
class LoocvPrimal;

template <typename T>
class LoocvDual;

template <typename T>
class FixLambda;

template <typename T>
class FixSigLam;

template <typename T>
class CalibrateSGD;

template <typename T>
class Siglam;

template <typename T>
class SiglamHo;

template <typename T>
class HoPrimal;

template <typename T>
class HoDual;

    /**
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
     * \return adds to opt the field paramsel 
     */

    virtual void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt) = 0;

    class BadParamSelectionCreation : public std::logic_error {
    public:
      BadParamSelectionCreation(std::string type)
      : logic_error("Cannot create type " + type) {}
    };
    static ParamSelection<T>*
    factory(const std::string& id) throw(BadParamSelectionCreation)
    {
      if(id == "loocvprimal")
        return new LoocvPrimal<T>;
      else if(id == "loocvdual")
        return new LoocvDual<T>;
      else if(id == "fixlambda")
        return new FixLambda<T>;
      else if(id == "calibratesgd")
        return new CalibrateSGD<T>;
      else if(id == "siglam")
        return new Siglam<T>;
      else if(id == "siglamho")
        return new SiglamHo<T>;
      else if(id == "hodual")
        return new HoDual<T>;
      else if(id == "hoprimal")
        return new HoPrimal<T>;
      else if(id == "fixsiglam")
        return new FixSigLam<T>;
      else
        throw BadParamSelectionCreation(id);
    }
};

}

#endif // _GURLS_PARAMSEL_H_
