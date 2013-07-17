#include "gurls++/gprwrapper.h"

#include "gurls++/gurls.h"
#include "gurls++/predkerneltraintest.h"
#include "gurls++/primal.h"
#include "gurls++/dual.h"
#include "gurls++/exceptions.h"

namespace gurls
{

template <typename T>
GPRWrapper<T>::GPRWrapper(const std::string &name): KernelWrapper<T>(name) { }

template <typename T>
void GPRWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    this->opt->removeOpt("split");
    this->opt->removeOpt("optimizer");

    OptTaskSequence *seq = new OptTaskSequence();
    GurlsOptionsList * process = new GurlsOptionsList("processes", false);
    OptProcess* process1 = new OptProcess();
    process->addOpt("one", process1);
    this->opt->addOpt("seq", seq);
    this->opt->addOpt("processes", process);
    this->opt->printAll();

    if(this->kType == KernelWrapper<T>::RBF) //could it be a different kernel?
    {
            *seq << "split:ho" << "paramsel:siglamhogpregr" << "kernel:rbf"; //is the order important...?
            *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
    }
    else return;

    *seq << "optimizer:rlsgpregr";
    *process1 << GURLS::computeNsave;

    GURLS G;
    G.run(X, y, *(this->opt), "one");

}

template <typename T>
gMat2D<T>* GPRWrapper<T>::eval(const gMat2D<T> &X)
{

    Prediction<T> *pred= new PredGPRegr<T>();
    PredKernelTrainTest<T> predkTrainTest;

    gMat2D<T> empty;

    this->opt->removeOpt("predkernel", true);
    this->opt->addOpt("predkernel", predkTrainTest.execute(X, empty, *(this->opt)));

    GurlsOptionsList* result = GurlsOptionsList::dynacast(pred->execute(X, empty, *(this->opt)));
    result->printAll();
    OptMatrix<gMat2D<T> > *means=result->getOptAs< OptMatrix<gMat2D<T> > >("means");
//     means->detachValue();
    gMat2D<T>* ret = &(means->getValue());
//     delete pred;
//     delete result;
//     delete means;
    return ret;
}

}
