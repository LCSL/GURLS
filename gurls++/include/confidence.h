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

#ifndef _GURLS_CONFIDENCE_H_
#define _GURLS_CONFIDENCE_H_


#include <stdexcept>
#include "optlist.h"

namespace gurls
{

template<typename T>
class ConfBoltzman;

template<typename T>
class ConfBoltzmanGap;

template<typename T>
class ConfGap;

template<typename T>
class ConfMaxScore;

/**
 * \ingroup Exceptions
 *
 * \brief BadConfidenceCreation is thrown when \ref factory tries to generate an unknown confidence method
 */
class BadConfidenceCreation : public std::logic_error
{
public:
    /**
     * Exception constructor.
     */
    BadConfidenceCreation(std::string type): logic_error("Cannot create type " + type) {}
};

/**
 * \ingroup Confidence
 * \brief Confidence is the class that computes a confidence score for the predicted labels
 */
template<typename T>
class Confidence
{
public:
    /**
     * Computes a confidence score for the predicted labels specified in the field pred of opt
     * \param X not used
     * \param Y not used
     * \param opt options with the following:
     *  - pred (settable with the class Prediction and its subclasses)
     *
     * \return ret, a GurlsOptionList with the following fields:
     *  - confidence = array containing the confidence score for each row of the field pred of opt.
     *  - labels = array containing predicted class for each row of the field pred of opt.
     */
    virtual GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) = 0;

    /**
     * Factory function returning a pointer to the newly created object.
     *
     * \warning The returned pointer is a plain, un-managed pointer. The calling
     * function is responsible of deallocating the object.
     */
    static Confidence<T> *factory(const std::string& id) throw(BadConfidenceCreation)
    {
        if(id == "boltzman")
            return new ConfBoltzman<T>;
        if(id == "boltzmangap")
            return new ConfBoltzmanGap<T>;
        if(id == "gap")
                return new ConfGap<T>;
        if(id == "maxscore")
            return new ConfMaxScore<T>;
        throw BadConfidenceCreation(id);
    }
};


}

#endif // _GURLS_CONFIDENCE_H_

