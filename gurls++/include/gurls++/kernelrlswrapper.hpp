
#include "gurls++/gurls.h"
#include "gurls++/kernelrlswrapper.h"
#include "gurls++/predkerneltraintest.h"
#include "gurls++/primal.h"
#include "gurls++/dual.h"
#include "gurls++/exceptions.h"

namespace gurls
{

template <typename T>
KernelRLSWrapper<T>::KernelRLSWrapper(const std::string &name): KernelWrapper<T>(name) { }

template <typename T>
void KernelRLSWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    this->opt->removeOpt("split");
    this->opt->removeOpt("optimizer");
    this->opt->removeOpt("seq");
    this->opt->removeOpt("processes");
    this->opt->removeOpt("todisk");

    const unsigned long nlambda = static_cast<unsigned long>(this->opt->getOptAsNumber("nlambda"));
    const unsigned long nsigma = static_cast<unsigned long>(this->opt->getOptAsNumber("nsigma"));
	
	unsigned int n=X.rows();
	unsigned int d=X.cols();
    OptTaskSequence *seq = new OptTaskSequence();
    GurlsOptionsList * process = new GurlsOptionsList("processes", false);
    OptProcess* process1 = new OptProcess();
    process->addOpt("one", process1);
    this->opt->addOpt("seq", seq);
    this->opt->addOpt("processes", process);
    this->opt->addOpt("todisk", new OptNumber(0));

    if(this->kType == KernelWrapper<T>::LINEAR)
    {
        if(nlambda > 1ul)
        {
			if(n>d){
				*seq << "split:ho" << "paramsel:hoprimal";
				*process1 << GURLS::computeNsave << GURLS::computeNsave;}
			else{
				*seq << "split:ho" << "kernel:linear" << "paramsel:hodual";
				*process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
			}
        }
        else if(nlambda == 1ul)
        {
            if(this->opt->hasOpt("paramsel.lambdas"))
            { if(n<=d){
                *seq << "kernel:linear";
				*process1 << GURLS::computeNsave;}
             }
            else
                throw gException("Please set a valid value for the regularization parameter, calling setParam(value)");
        }
        else
            throw gException("Please set a valid value for NParam, calling setNParam(value)");

    }
    else if(this->kType == KernelWrapper<T>::RBF)
    {
        if(nlambda > 1ul)
        {
            if(nsigma > 1ul)
            {
                *seq << "split:ho" << "paramsel:siglamho" << "kernel:rbf";
                *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
            }
            else if(nsigma == 1ul)
            {
                if(this->opt->hasOpt("paramsel.sigma"))
                {
                    *seq << "split:ho" << "kernel:rbf" << "paramsel:hodual";
                    *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
                }
                else
                    throw gException("Please set a valid value for the kernel parameter, calling setSigma(value)");
            }
            else
                throw gException("Please set a valid value for NSigma, calling setNSigma(value)");
        }
        else if(nlambda == 1ul)
        {
            if(nsigma == 1ul)
            {
                if(this->opt->hasOpt("paramsel.sigma") && this->opt->hasOpt("paramsel.lambdas"))
                {
                    *seq << "kernel:rbf";
                    *process1 << GURLS::computeNsave;
                }
                else
                    throw gException("Please set a valid value for kernel and regularization parameters, calling setParam(value) and setSigma(value)");
            }
            else
                throw gException("Please set a valid value for NSigma, calling setNSigma(value)");
        }
        else
            throw gException("Please set a valid value for NParam, calling setNParam(value)");
    }

	if (this->kType == KernelWrapper<T>::LINEAR && n>d)
		*seq << "optimizer:rlsprimal";
	else
		*seq << "optimizer:rlsdual";

    *process1 << GURLS::computeNsave;

    GURLS G;
    G.run(X, y, *(this->opt), "one");
	
}

template <typename T>
gMat2D<T>* KernelRLSWrapper<T>::eval(const gMat2D<T> &X)
{
	
    if(!this->trainedModel())
        throw gException("Error, Train Model First");

    Prediction<T> *pred;
    PredKernelTrainTest<T> predkTrainTest;

    gMat2D<T> empty;
	
    switch (this->kType)
    {
    case KernelWrapper<T>::LINEAR:
        pred = new PredPrimal<T>();
        break;
    case KernelWrapper<T>::RBF:
        pred = new PredDual<T>();
        this->opt->removeOpt("predkernel");
        this->opt->addOpt("predkernel", predkTrainTest.execute(X, empty, *(this->opt)));
		break;
	default:
		throw gException("Kernel not recognized");
    }
    OptMatrix<gMat2D<T> >* result = OptMatrix<gMat2D<T> >::dynacast(pred->execute(X, empty, *(this->opt)));
    delete pred;
    result->detachValue();

    gMat2D<T>* ret = &(result->getValue());
    delete result;
    return ret;
}

}
