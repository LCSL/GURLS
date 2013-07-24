#include "gurls++/gprwrapper.h"

#include "gurls++/gurls.h"
#include "gurls++/predkerneltraintest.h"
#include "gurls++/dual.h"
#include "gurls++/exceptions.h"
#include "gurls++/normzscore.h"
#include "gurls++/normtestzscore.h"

namespace gurls
{

template <typename T>
GPRWrapper<T>::GPRWrapper(const std::string &name): KernelWrapper<T>(name), norm(NULL) { }

template <typename T>
GPRWrapper<T>::~GPRWrapper()
{
    if(norm != NULL)
        delete norm;
}

template <typename T>
void GPRWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    this->opt->removeOpt("split");
    this->opt->removeOpt("optimizer");

    if(norm != NULL)
    {
        delete norm;
        norm = NULL;
    }

    NormZScore<T> taskNorm;

    norm = taskNorm.execute(X, y, *(this->opt));



    const unsigned long nlambda = static_cast<unsigned long>(this->opt->getOptAsNumber("nlambda"));
    const unsigned long nsigma = static_cast<unsigned long>(this->opt->getOptAsNumber("nsigma"));


    OptTaskSequence *seq = new OptTaskSequence();
    GurlsOptionsList * process = new GurlsOptionsList("processes", false);
    OptProcess* process1 = new OptProcess();
    process->addOpt("one", process1);
    this->opt->addOpt("seq", seq);
    this->opt->addOpt("processes", process);
//    this->opt->printAll();

    if(this->kType == KernelWrapper<T>::LINEAR)
    {
        if(nlambda > 1ul)
        {
            *seq << "split:ho" << "kernel:linear" << "paramsel:hogpregr";
            *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
        }
        else if(nlambda == 1ul)
        {
            if(this->opt->hasOpt("paramsel.lambdas"))
            {
                *seq << "kernel:linear";
                *process1 << GURLS::computeNsave;
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
                *seq << "split:ho" << "paramsel:siglamhogpregr" << "kernel:rbf";
                *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
            }
            else if(nsigma == 1ul)
            {
                if(this->opt->hasOpt("paramsel.sigma"))
                {
                    *seq << "split:ho" << "kernel:rbf" << "paramsel:hogpregr";
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

    *seq << "optimizer:rlsgpregr";
    *process1 << GURLS::computeNsave;

    GURLS G;
    G.run(norm->getOptValue<OptMatrix<gMat2D<T> > >("X"), norm->getOptValue<OptMatrix<gMat2D<T> > >("Y"), *(this->opt), "one");

}

template <typename T>
gMat2D<T>* GPRWrapper<T>::eval(const gMat2D<T> &X)
{
    gMat2D<T> vars;

    return eval(X, vars);
}

template <typename T>
gMat2D<T>* GPRWrapper<T>::eval(const gMat2D<T> &X, gMat2D<T> &vars)
{
    gMat2D<T> empty;

    NormTestZScore<T> normTask;
    PredGPRegr<T> predTask;


    GurlsOptionsList *normX = normTask.execute(X, empty, *(norm));

    gMat2D<T> &Xresc = normX->getOptValue<OptMatrix<gMat2D<T> > >("X");

    PredKernelTrainTest<T> predkTrainTest;
    this->opt->removeOpt("predkernel");
    this->opt->addOpt("predkernel", predkTrainTest.execute(Xresc, empty, *(this->opt)));

    GurlsOptionsList *pred = predTask.execute(Xresc, empty, *(this->opt));

    delete normX;

    OptMatrix<gMat2D<T> >* pmeans = pred->getOptAs<OptMatrix<gMat2D<T> > >("means");
    pmeans->detachValue();

    gMat2D<T> &predMeans = pmeans->getValue();
    gMat2D<T> &predVars = pred->getOptValue<OptMatrix<gMat2D<T> > >("vars");

    const unsigned long n = predMeans.rows();
    const unsigned long t = predMeans.cols();

    T* column = predMeans.getData();
    const T* std_it = norm->getOptValue<OptMatrix<gMat2D<T> > >("stdY").getData();
    const T* mean_it = norm->getOptValue<OptMatrix<gMat2D<T> > >("meanY").getData();
    const T* pvars_it = predVars.getData();

    vars.resize(n, t);

    T* vars_it = vars.getData();

    for(unsigned long i=0; i<t; ++i, column+=n, ++std_it, ++mean_it, vars_it+=n)
    {
        scal(n, *std_it, column, 1);
        axpy(n, (T)1.0, mean_it, 0, column, 1);

        copy(vars_it, pvars_it, n);
        scal(n, (*std_it)*(*std_it), vars_it, 1);
    }

    delete pred;

    return &predMeans;
}

}
