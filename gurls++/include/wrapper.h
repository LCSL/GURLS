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
public:

    enum ProblemType{CLASSIFICATION, REGRESSION};


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
    virtual const GurlsOptionsList& getOpt() const;

    virtual void saveModel(const std::string &fileName);
    virtual void loadModel(const std::string &fileName);

//    virtual void exportModel(const std::string &fileName);
//    virtual void importModel(const std::string &fileName);

    virtual void setNparams(unsigned long value);
    virtual void setParam(double value);

    virtual void setSplitProportion(double value);
    virtual void setProblemType(ProblemType value);


protected:
    /**
      * Checks if model has already been trained
      */
    virtual bool trainedModel();

    std::string name;       ///< Name of the options structure
    GurlsOptionsList *opt;  ///< Options structure where information about initial training is stored

    ProblemType probType;

};

template<typename T>
class KernelWrapper : public GurlsWrapper<T>
{
public:
    enum KernelType{RBF, LINEAR, CHISQUARED};

    KernelWrapper(const std::string& name);

    virtual void train(const gMat2D<T> &X, const gMat2D<T> &y) = 0;

    virtual void setKernelType(KernelType value);
    virtual void setSigma(double value);
};

}

#include"wrapper.hpp"

#endif //GURLS_WRAPPER_H
