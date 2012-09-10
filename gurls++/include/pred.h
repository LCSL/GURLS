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


#ifndef _GURLS_PRED_H_
#define _GURLS_PRED_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <exception>
#include <stdexcept>

#include "gmath.h"
#include "options.h"
#include "optlist.h"


namespace gurls {

template <typename T>
class PredPrimal;

template <typename T>
class PredDual;

template <typename T>
class PredGPRegr;

/**
 * \ingroup Prediction
 * \brief Prediction is the class that computes predictions
 */

template <typename T>
class Prediction
{
public:
    /**
     * Computes predictions of the classifier stored in the field rls of opt on the samples passed in the X matrix.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the different required fields based on the sub-class
     *
     * \return adds the field pred to opt
     */
    virtual GurlsOption* execute( const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) = 0;

    /**
     * \ingroup Exceptions
     *
     * \brief BadPredictionCreation is thrown when \ref factory tries to generate an unknown prediction method
     */
    class BadPredictionCreation : public std::logic_error
    {
    public:
        /**
         * Exception constructor.
         */
        BadPredictionCreation(std::string type)
            : logic_error("Cannot create type " + type) {}
    };

    /**
     * Factory function returning a pointer to the newly created object.
     *
     * \warning The returned pointer is a plain, un-managed pointer. The calling
     * function is responsible of deallocating the object.
     */
    static Prediction<T>* factory(const std::string& id) throw(BadPredictionCreation)
    {
        if(id == "primal")
            return new PredPrimal<T>;
        else if(id == "dual")
            return new PredDual<T>;
        else if(id == "gpregr")
            return new PredGPRegr<T>;
        else
            throw BadPredictionCreation(id);
    }

};

}
#endif // _GURLS_PRED_H_
