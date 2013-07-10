#include "rlswrapper.h"

#include "splitho.h"
#include "hoprimal.h"
#include "rlsprimal.h"
#include "primal.h"

namespace gurls
{

template <typename T>
RLSWrapper<T>::RLSWrapper(const std::string &name): GurlsWrapper<T>(name) { }

template <typename T>
void RLSWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    this->opt->removeOpt("split");
    this->opt->removeOpt("optimizer");

    const unsigned long nlambda = static_cast<unsigned long>(this->opt->getOptAsNumber("nlambda"));
    if(nlambda != 1)
    {
        SplitHo<T> splitTask;
        this->opt->addOpt("split", splitTask.execute(X, y, *(this->opt)));

        ParamSelHoPrimal<T> paramselTask;
        GurlsOption * result = paramselTask.execute(X, y, *(this->opt));
        this->opt->removeOpt("paramsel");
        this->opt->addOpt("paramsel", result);
    }
    else
    {
        if(!this->opt->hasOpt("paramsel.lambdas"))
            throw gException("Please set a valid value for the regularization parameter, calling setParam(value)");
    }

    RLSPrimal<T> optimizerTask;
    this->opt->addOpt("optimizer", optimizerTask.execute(X, y, *(this->opt)));
}

template <typename T>
gMat2D<T>* RLSWrapper<T>::eval(const gMat2D<T> &X)
{
    PredPrimal<T> predTask;
    gMat2D<T> empty;
    OptMatrix<gMat2D<T> >* result = predTask.execute(X, empty, *(this->opt));

    gMat2D<T>* pred = &(result->getValue());

    result->detachValue();
    delete result;

    return pred;
}

}
