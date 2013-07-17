/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
  * Authors: Matteo Santoro, Elena Ceseracciu
  TODO: add copy policy
 */

#ifndef GURLS_GPRWRAPPER_H
#define GURLS_GPRWRAPPER_H

#include "gurls++/wrapper.h"

namespace gurls
{

/**
  * \ingroup Wrappers
  * \brief GPRWrapper is the sub-class of GurlsWrapper that implements ...
  */
template<typename T>
class GPRWrapper: public KernelWrapper<T>
{
public:
    /**
      * Constructor
      *
      * \param name Name of the option's structure that will be initialized
      */
    GPRWrapper(const std::string& name);

    /**
      * Initial parameter selection and training
      *
      * \param X Input data matrix
      * \param Y Labels matrix
      */
    void train(const gMat2D<T> &X, const gMat2D<T> &y);

    /**
      * Estimates label for an input matrix
      *
      * \param[in] X Input matrix
      * \returns Matrix of predicted labels
      */
    gMat2D<T>* eval(const gMat2D<T> &X);
};

}

#include "gprwrapper.hpp"

#endif //GURLS_KERNELRLSWRAPPER_H
