#ifndef GURLS_WRAPPER_H
#define GURLS_WRAPPER_H

#include "gvec.h"
#include "gmat2d.h"
#include "optlist.h"

namespace gurls
{

class GurlsWrapper
{
protected:
    typedef double T;

public:
    GurlsWrapper(const std::string& name);
    ~GurlsWrapper();

    virtual void train(const gMat2D<T> &X, const gMat2D<T> &y) = 0;
    virtual void update(const gVec<T> &X, const gVec<T> &y) = 0;
    virtual T eval(const gVec<T> &X, unsigned long *index = NULL);
    virtual gMat2D<T>* eval(const gMat2D<T> &X) = 0;

    const GurlsOptionsList& getOpt() const;

protected:
    virtual bool trainedModel();

    std::string name;
    GurlsOptionsList *opt;
};


class RecursiveRLSWrapper: public GurlsWrapper
{
public:
    RecursiveRLSWrapper(const std::string& name);

    void train(const gMat2D<T> &X, const gMat2D<T> &y);
    void update(const gVec<T> &X, const gVec<T> &y);
    gMat2D<T>* eval(const gMat2D<T> &X);
    using GurlsWrapper::eval;

};


class RecursiveRLSRetrainWrapper: public RecursiveRLSWrapper
{
public:
    RecursiveRLSRetrainWrapper(const std::string& name);

    void train(const gMat2D<T> &X, const gMat2D<T> &y);
    void update(const gVec<T> &X, const gVec<T> &y);
    void retrain();

protected:
    unsigned long nTot;
};

}
#endif //GURLS_WRAPPER_H
