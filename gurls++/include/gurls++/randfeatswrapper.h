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

#ifndef GURLS_RANDFEATSWRAPPER_H
#define GURLS_RANDFEATSWRAPPER_H

#include "gurls++/rlswrapper.h"

namespace gurls
{

/**
  * \ingroup Wrappers
  * \brief RLSWrapper is the sub-class of GurlsWrapper that implements Regularized Least Squares with a linear model.
  * The regularization parameter is estimated by the wrapper (default) or explicitely given in input by the user via method train().
  * The eval() method estimates the output for new data.
  *
  */
template<typename T>
class RandomFeaturesWrapper: public RLSWrapper<T>
{
public:
    /**
      * Constructor
      *
      * \param name Name of the option's structure that will be initialized
      */
    RandomFeaturesWrapper(const std::string& name);

    /**
      * Destructor
      */
    ~RandomFeaturesWrapper();

    /**
      * Initial parameter selection and training
      *
      * \param[in] X Input data matrix
      * \param[in] Y Labels matrix
      */
    void train(const gMat2D<T> &X, const gMat2D<T> &y);


    /**
      * Estimates label for an input matrix
      *
      * \param[in] X Input matrix
      * \returns Matrix of predicted labels
      */
    gMat2D<T> *eval(const gMat2D<T> &X);

    /**
      *
      * \param value
      */
    void setNRandFeats(unsigned long value);

protected:
    gMat2D<T> *W;
};

}

#include "randfeatswrapper.hpp"

#endif //GURLS_RANDFEATSWRAPPER_H
