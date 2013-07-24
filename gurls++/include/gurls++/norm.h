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

#ifndef _GURLS_NORM_H_
#define _GURLS_NORM_H_

#include <stdexcept>
#include <cstring>
#include <cmath>

#include "gurls++/gvec.h"
#include "gurls++/gmat2d.h"
#include "gurls++/exceptions.h"

#include "gurls++/optlist.h"

namespace gurls
{

static const std::string L0norm = "l0";     ///< String identifying l0-norm
static const std::string L1norm = "l1";     ///< String identifying l1-norm
static const std::string L2norm = "l2";     ///< String identifying l2-norm
static const std::string LInfnorm = "inf";  ///< String identifying infinity-norm

/**
  * Computes a vector norm specified by \a type parameter
  */
template<typename T>
T norm(const gVec<T>& x, std::string type = "l2"){
    T nrm = 0;
    const T* ptr = x.getData();

    if (type.compare(L0norm) == 0){
        for(unsigned long i =0;i < x.getSize(); ++i){
            if (*ptr++ != 0)
                nrm++;
        }
        return nrm;
    }else if (type.compare(L1norm) == 0) {
        for(unsigned long i =0; i < x.getSize(); ++i){
            nrm+= static_cast<T>(std::abs(*ptr++));
        }
        return nrm;
    }else if (type.compare(L2norm) == 0) {
        for(unsigned long i =0; i < x.getSize(); ++i, ++ptr){
            nrm+= std::pow(*ptr, 2);
        }
        return static_cast<T>(std::sqrt(nrm));
    }else if (type.compare(LInfnorm) == 0 ) {
        throw gException("Sorry. You required the computation of LInfnorm but this functionality has not been implemented yet!");
    }else {
        throw gException("Unknown norm type.");
    }

    return nrm;
}

/**
  * Computes a matrix norm specified by \a type parameter
  */
template<typename T>
T norm(const gMat2D<T>& A, std::string type = "l2"){

    throw gException("Sorry. You required to compute the norm of a matrix but this functionality has not been implemented yet!");
}

template<typename T>
class NormL2;

template<typename T>
class NormZScore;

template<typename T>
class NormTestZScore;

/**
 * \ingroup Exceptions
 *
 * \brief BadNormCreation is thrown when \ref factory tries to generate an unknown norm
 */
class BadNormCreation : public std::logic_error
{
public:

    /**
     * Exception constructor.
     */
    BadNormCreation(const std::string& type): logic_error("Cannot create type " + type) {}
};

/**
 * \ingroup Norms
 * \brief Norm is a class that spherifies the data
 */
template<typename T>
class Norm
{
public:

    /**
     * Spherifies the data
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the different required fields based on the sub-class
     * \return spherified input data matrix
     */
    virtual GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) = 0;

    /**
     * Factory function returning a pointer to the newly created object.
     *
     * \warning The returned pointer is a plain, un-managed pointer. The calling
     * function is responsible of deallocating the object.
     */
    static Norm<T> *factory(const std::string& id) throw(BadNormCreation)
    {
        if(id == "l2")
            return new NormL2<T>;
        if(id == "zscore")
            return new NormZScore<T>;
        if(id == "testzscore")
            return new NormTestZScore<T>;

        throw BadNormCreation(id);
    }
};

}

#endif // _GURLS_NORM_H_

