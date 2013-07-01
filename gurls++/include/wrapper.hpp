
#include "gurls.h"
#include "optlist.h"
#include "wrapper.h"

namespace gurls
{

template <typename T>
GurlsWrapper<T>::GurlsWrapper(const std::string& name):opt(NULL), name(name)
{
    opt = new GurlsOptionsList(name, true);

    setSplitProportion(0.2);
    setNparams(20);
    setProblemType(CLASSIFICATION);
}

template <typename T>
GurlsWrapper<T>::~GurlsWrapper()
{
    delete opt;
}

template <typename T>
T GurlsWrapper<T>::eval(const gVec<T> &X, unsigned long *index)
{
    if(!trainedModel())
        throw gException("Error, Train Model First");

    gMat2D<T>X_mat(1, X.getSize());
    copy(X_mat.getData(), X.getData(), X.getSize());

    gMat2D<T>* pred_mat = eval(X_mat);

    const T* pred = pred_mat->getData();
    const unsigned long size = pred_mat->getSize();

    const T* max = std::max_element(pred, pred+size);
    T ret = *max;
    if(index != NULL)
        *index = max-pred;

    delete pred_mat;
    return ret;
}

template <typename T>
const GurlsOptionsList &GurlsWrapper<T>::getOpt() const
{
    return *opt;
}

template <typename T>
void GurlsWrapper<T>::saveModel(const std::string &fileName)
{
    opt->save(fileName);
}

template <typename T>
void GurlsWrapper<T>::loadModel(const std::string &fileName)
{
    opt->load(fileName);
}

template <typename T>
void GurlsWrapper<T>::setNparams(unsigned long value)
{
    opt->getOptValue<OptNumber>("nlambda") = value;

    if(opt->hasOpt("paramsel.lambdas") && value > 1.0)
    {
        std::cout << "Warning: ignoring previous values of the regularization parameter" << std::endl;
        opt->getOptAs<GurlsOptionsList>("paramsel")->removeOpt("lambdas");
    }
}

template <typename T>
void GurlsWrapper<T>::setParam(double value)
{
    if(!opt->hasOpt("paramsel"))
        opt->addOpt("paramsel", new GurlsOptionsList("paramsel"));

    if(opt->hasOpt("paramsel.lambdas"))
        opt->getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas").getData()[0] = (T)value;
    else
    {
        gMat2D<T> * lambdas = new gMat2D<T>(1,1);
        lambdas->getData()[0] = (T)value;
        opt->getOptAs<GurlsOptionsList>("paramsel")->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*lambdas));
    }

    setNparams(1);
}

template <typename T>
void GurlsWrapper<T>::setSplitProportion(double value)
{
    opt->getOptValue<OptNumber>("hoproportion") = value;
}

template <typename T>
void GurlsWrapper<T>::setProblemType(typename GurlsWrapper::ProblemType value)
{
    probType = value;

    opt->getOptValue<OptString>("hoperf") = (value == CLASSIFICATION)? "macroavg": "rmse";
}

template <typename T>
bool GurlsWrapper<T>::trainedModel()
{
    return opt->hasOpt("optimizer");
}


template <typename T>
KernelWrapper<T>::KernelWrapper(const std::string &name): GurlsWrapper<T>(name), kType(RBF)
{
//    GurlsOptionsList *kernel = new GurlsOptionsList("kernel");
//    this->opt->addOpt("kernel", kernel);
//    kernel->addOpt("type", "rbf");

    GurlsOptionsList *paramsel = new GurlsOptionsList("paramsel");
    this->opt->addOpt("paramsel", paramsel);
}

template <typename T>
void KernelWrapper<T>::setKernelType(typename KernelWrapper::KernelType value)
{
    kType = value;

//    std::string &type = this->opt->template getOptValue<OptString>("kernel.type");

//    switch(value)
//    {
//    case RBF:
//        type = std::string("rbf");
//    case LINEAR:
//        type = std::string("linear");
//    case CHISQUARED:
//        type = std::string("chisquared");
//    }
}

template <typename T>
void KernelWrapper<T>::setSigma(double value)
{
    if(this->opt->hasOpt("paramsel.sigma"))
        this->opt->template getOptValue<OptNumber>("paramsel.sigma") = value;
    else
    {
        GurlsOptionsList* paramsel = this->opt->template getOptAs<GurlsOptionsList>("paramsel");
        paramsel->addOpt("sigma", new OptNumber(value));
    }

    setNSigma(1);
}

template <typename T>
void KernelWrapper<T>::setNSigma(unsigned long value)
{
    this->opt->template getOptValue<OptNumber>("nsigma") = value;

    if(this->opt->hasOpt("paramsel.sigma") && value > 1.0)
    {
        std::cout << "Warning: ignoring previous values of the kernel parameter" << std::endl;
        this->opt->template getOptAs<GurlsOptionsList>("paramsel")->removeOpt("sigma");
    }
}

}
