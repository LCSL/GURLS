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

#include "gurls++/gvec.h"
#include "gurls++/gmat2d.h"
#include "gurls++/optlist.h"

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
      * \param[in] X Input data matrix
      * \param[in] Y Labels matrix
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
      * Estimates performance for an input matrix
      *
      * \param[in] X Input matrix
      * \returns Matrix of predicted labels
      */
    virtual gMat2D<T>* perf(const gMat2D<T> &y, gMat2D<T> &pred, const std::string perfname);

    /**
      * Returns a const reference to the options structure
      */
    virtual const GurlsOptionsList& getOpt() const;

    /**
      * Saves the computed model to file
      *
      * \param fileName name of the file where data will be saved
      */
    virtual void saveModel(const std::string &fileName);

	 /**
      * Sets the file in which the model will be saved
      *
      * \param fileName name of the file where data will be saved
      */
	virtual void setSavefile(const std::string &fileName);

    /**
      * Loads a computed model from a file
      *
      * \param fileName name of the file containing the data to load
      */
    virtual void loadModel(const std::string &fileName);

	virtual void loadOpt(GurlsOptionsList &opt);

//    virtual void exportModel(const std::string &fileName);
//    virtual void importModel(const std::string &fileName);

    /**
      *
      * \param[in] value
      */
    virtual void setNparams(unsigned long value);
    /**
      *
      * \param[in] value
      */
    virtual void setParam(double value);	
    /**
      *
      * \param[in] value
      */
    virtual void setSplitProportion(double value);
    /**
      *
      * \param[in] value
      */
    virtual void setProblemType(ProblemType value);
    /**
      *
      * \param[out] problem type
      */
	ProblemType getProblemType();
    /**
      * Estimates problem type
      *
      * \param[in] X Input data matrix
      * \param[in] Y Labels matrix      
	  * \param[out] probType estimated problem type
      * \returns Estimated problem type
      */	
	ProblemType problemTypeFromData( const gMat2D<T> &X, const gMat2D<T> &y);

protected:
    /**
      * Checks if model has already been trained
      */
    virtual bool trainedModel();

    std::string name;       ///< Name of the options structure
    GurlsOptionsList *opt;  ///< Options structure where information about initial training is stored

    ProblemType probType;   ///< Problem type

};

/**
  * \ingroup Wrappers
  * \brief KernelWrapper is the base class for all gurls++ wrappers
  */
template<typename T>
class KernelWrapper : public GurlsWrapper<T>
{
public:
    enum KernelType{RBF, LINEAR/*, CHISQUARED*/};                   ///<

    /**
      * Constructor
      *
      * \param name Name of the options structure that will be initialized
      */
    KernelWrapper(const std::string& name);

    /**
      * Initial parameter selection and training
      *
      * \param X Input data matrix
      * \param Y Labels matrix
      */
    virtual void train(const gMat2D<T> &X, const gMat2D<T> &y) = 0;

    /**
      *
      * \param value
      */
    virtual void setKernelType(KernelType value);

    /**
      *
      * \param value
      */
    virtual void setSigma(double value);

    /**
      *
      * \param value
      */
    virtual void setNSigma(unsigned long value);

protected:
    KernelType kType;   ///< Kernel type used in train and eval
};

}

#include"wrapper.hpp"

#endif //GURLS_WRAPPER_H
