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

#ifndef GURLS_WRAPPER_H
#define GURLS_WRAPPER_H

#include "gvec.h"
#include "gmat2d.h"
#include "optlist.h"

namespace gurls
{

/**
  * \ingroup Wrappers
  * \brief GurlsWrapper is the base class for all gurls++ wrappers
  */
template<typename T>
class GurlsWrapper
{
protected:
//    typedef double T;   ///< Data type for matrices cells

public:
    /**
      * Constructor
      *
      * \param name Name of the options structure that will be initialized
      */
    GurlsWrapper(const std::string& name);

    /**
      * Destructor
      */
    ~GurlsWrapper();

    /**
      * Initial parameter selection and training
      *
      * \param X Input data matrix
      * \param Y Labels matrix
      */
    virtual void train(const gMat2D<T> &X, const gMat2D<T> &y) = 0;

    /**
      * Estimator update
      *
      * \param X Input data vector
      * \param Y Labels vector
      */
    virtual void update(const gVec<T> &X, const gVec<T> &y) = 0;

    /**
      * Estimates label for a new input point
      *
      * \param[in] X Input point
      * \param[out] index Index of the estimated label
      * \returns Estimated label
      */
    virtual T eval(const gVec<T> &X, unsigned long *index = NULL);

    /**
      * Estimates label for an input matrix
      *
      * \param[in] X Input matrix
      * \returns Matrix of predicted labels
      */
    virtual gMat2D<T>* eval(const gMat2D<T> &X) = 0;

    /**
      * Returns a const reference to the options structure
      */
    const GurlsOptionsList& getOpt() const;

protected:
    /**
      * Checks if model has already been trained
      */
    virtual bool trainedModel();

    std::string name;       ///< Name of the options structure
    GurlsOptionsList *opt;  ///< Options structure where information about initial training is stored

};

/**
  * \ingroup Wrappers
  * \brief RecursiveRLSWrapper is the sub-class of GurlsWrapper that implements recursive update
  * of the RLS estimator without retraining.
  *
  * Initial parameter selection and training are carried out on a initial set of samples. The
  * computation of the RLS estimator is carried out by the class RLSPrimalRecInit,
  * which stores all information necessary for efficient recursive update in the options structure.
  * Once the information about initial training is stored in the options structure, given a
  * new input–output pair, the RLS estimator can be efficiently updated via the method update().
  * Every time a new input-output pair is available, method update() can be invoked again.
  * Finally, the eval() method can be used on test data.
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
      * \param X Input data matrix
      * \param Y Labels matrix
      */
    void train(const gMat2D<T> &X, const gMat2D<T> &y);

    /**
      * Estimator update
      *
      * \param X Input data vector
      * \param Y Labels vector
      */
    void update(const gVec<T> &X, const gVec<T> &y);

    /**
      * Estimates label for an input matrix
      *
      * \param[in] X Input matrix
      * \returns Matrix of predicted labels
      */
    gMat2D<T>* eval(const gMat2D<T> &X);


    using GurlsWrapper<T>::eval;

};

/**
  * \ingroup Wrappers
  * \brief RecursiveRLSRetrainWrapper is the sub-class of GurlsWrapper that implements recursive update
  * of the RLS estimator with retraining capability.
  *
  * Initial parameter selection and training are carried out on a initial set of samples. The
  * computation of the RLS estimator is carried out by the class RLSPrimalRecInit,
  * which stores all information necessary for efficient recursive update in the options structure.
  * Once the information about initial training is stored in the options structure, given a
  * new input–output pair, the RLS estimator can be efficiently updated via the method update().
  * Every time a new input-output pair is available, method update() can be invoked again. Parameter selection
  * and RLS estimation ( method retrain()) can be repeated after any number of online updates.
  * Finally, the eval() method can be used on test data.
  */
template<typename T>
class RecursiveRLSRetrainWrapper: public RecursiveRLSWrapper<T>
{
public:
    /**
      * Constructor
      *
      * \param name Name of the option's structure that will be initialized
      */
    RecursiveRLSRetrainWrapper(const std::string& name);

    /**
      * Initial parameter selection and training
      *
      * \param X Input data matrix
      * \param Y Labels matrix
      */
    void train(const gMat2D<T> &X, const gMat2D<T> &y);

    /**
      * Estimator update
      *
      * \param X Input data vector
      * \param Y Labels vector
      */
    void update(const gVec<T> &X, const gVec<T> &y);

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

#include"wrapper.hpp"

#endif //GURLS_WRAPPER_H
