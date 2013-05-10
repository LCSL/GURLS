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
class GurlsWrapper
{
protected:
    typedef double T;   ///< Data type for matrices cells

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
  * new set of input–output pairs (possibly also just one pair), the RLS estimator can be efficiently
  * updated via the method update(). Every time a new set of input-output
  * pairs is available, method update() can be invoked again.
  * Finally, the eval() method can be used on test data.
  */
class RecursiveRLSWrapper: public GurlsWrapper
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


    using GurlsWrapper::eval;

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
  * new set of input–output pairs (possibly also just one pair), the RLS estimator can be efficiently
  * updated via the method update(). Every time a new set of input-output
  * pairs is available, method update() can be invoked again. Parameter selection
  * and RLS estimation ( method retrain()) can be repeated after any number of online updates.
  * Finally, the eval() method can be used on test data.
  */
class RecursiveRLSRetrainWrapper: public RecursiveRLSWrapper
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
#endif //GURLS_WRAPPER_H
