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

#ifndef GURLS_RECRLSWRAPPER_H
#define GURLS_RECRLSWRAPPER_H

#include "wrapper.h"

namespace gurls
{

/**
  * \ingroup Wrappers
  * \brief RecursiveRLSWrapper is the sub-class of GurlsWrapper that implements recursive update
  * of the RLS estimator with retraining capability.
  *
  * Initial parameter selection and training are carried out on a initial set of samples via method train() 
  * which stores all information necessary for efficient recursive update in the options structure.
  * Once the information about initial training is stored, given a new input-output pair, 
  * the RLS estimator can be efficiently updated via the method update().
  * Every time a new input-output pair is available, method update() can be invoked again. Parameter selection and RLS estimation ( method retrain()) can be repeated after any number of online updates.
  * Finally, the eval() method estimates the output for new data.
  */
template<typename T>
class RecursiveRLSWrapper: public GurlsWrapper<T>
{
public:
    /**
      * Constructor
      *
      * \param name Name of the option's structure that will be initialized
      */
    RecursiveRLSWrapper(const std::string& name);

    /**
      * Initial parameter selection and training
      *
      * \param[in] X Input data matrix
      * \param[in] Y Labels matrix
      */
    void train(const gMat2D<T> &X, const gMat2D<T> &y);

    /**
      * Estimator update
      *
      * \param[in] X Input data vector
      * \param[in] Y Labels vector
      */
    void update(const gVec<T> &X, const gVec<T> &y);

    /**
      * Estimates label for an input matrix
      *
      * \param[in] X Input matrix
      * \returns Matrix of predicted labels
      */
    gMat2D<T>* eval(const gMat2D<T> &X);

    /**
      * Selection of the new regularization parameter.
      * \brief Selection is performed via hold-out validation using the subset
      * Xva,yva of the total training set as validation set.
      */
    void retrain();

protected:
    unsigned long nTot; ///< Total number of samples used for training
};

}

#include"recrlswrapper.hpp"

#endif //GURLS_RECRLSWRAPPER_H
