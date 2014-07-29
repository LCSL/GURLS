
#include "gurls++/gurls.h"
#include "gurls++/optlist.h"
#include "gurls++/wrapper.h"

namespace gurls
{

template <typename T>
GurlsWrapper<T>::GurlsWrapper(const std::string& name):opt(NULL), name(name)
{
    opt = new GurlsOptionsList(name, true);
	this->isowner=true;
    setSplitProportion(0.2);
    setNparams(20);
    setProblemType(CLASSIFICATION);
}

template <typename T>
GurlsWrapper<T>::~GurlsWrapper()
{
if(this->isowner)
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
void GurlsWrapper<T>::setSavefile(const std::string &fileName)
{
		opt->removeOpt("savefile");
		opt->addOpt("savefile",fileName);
}

template <typename T>
void GurlsWrapper<T>::loadOpt(const std::string &fileName)
{
	if(this->isowner)
		delete opt;
    opt = new GurlsOptionsList(name,false);
	this->isowner=true;
    opt->load(fileName);
}

template <typename T>
void GurlsWrapper<T>::loadOpt(GurlsOptionsList &optnew, bool owner)
{
	if(this->isowner)
		delete opt;
	if(owner)
	{
		opt = new GurlsOptionsList(optnew);
		this->isowner=true;
	}
	else
	{
		opt = &optnew;
		this->isowner=false;
	}
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
typename GurlsWrapper<T>::ProblemType GurlsWrapper<T>::getProblemType()
{
    return probType;
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

template <typename T>
typename GurlsWrapper<T>::ProblemType GurlsWrapper<T>::problemTypeFromData( const gMat2D<T> &X, const gMat2D<T> &y)
{
	// if 
	// max(max(abs(double(int32(y))-y))) > eps
	// probType = REGRESSION
	// else
	// probType = CLASSIFICATION

	for(unsigned long i=0; i<y.cols(); ++i)
		for (unsigned long j=0; j<y.rows(); ++j)
			if((y[j][i]-std::floor(y[j][i]))>0)
			{return REGRESSION;}
	return CLASSIFICATION;

}

template <typename T>
gMat2D<T>* GurlsWrapper<T>::perf(const gMat2D<T> &y, gMat2D<T> &pred, const std::string perfstring)
{
	std::string name=perfstring;
	if (perfstring=="macroavg")
		name=std::string("acc");
    Performance<T> *perf = Performance<T>::factory(perfstring);
    gMat2D<T> empty;

    this->opt->removeOpt("pred");
    this->opt->addOpt("pred", new OptMatrix<gMat2D<T> >(pred));

	GurlsOptionsList* perfList =perf->execute(empty, y, *(this->opt));
	gMat2D<T>* ret = &(perfList->getOptValue<OptMatrix<gMat2D<T> > >(name));

	delete perf;

	
    return ret;
}

}
