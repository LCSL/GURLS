
#include "gurls.h"
#include "optlist.h"
#include "wrapper.h"

namespace gurls
{
template <typename T>
GurlsWrapper<T>::GurlsWrapper(const std::string& name):opt(NULL), name(name) {}

template <typename T>
GurlsWrapper<T>::~GurlsWrapper()
{
    if(opt != NULL)
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
bool GurlsWrapper<T>::trainedModel()
{
    if(opt == NULL)
        return false;

    return opt->hasOpt("optimizer");
}


template <typename T>
RecursiveRLSWrapper<T>::RecursiveRLSWrapper(const std::string& name):GurlsWrapper<T>(name) {}

template <typename T>
void RecursiveRLSWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    if(this->opt != NULL)
        delete this->opt;

    this->opt = new GurlsOptionsList(this->name, true);

    OptTaskSequence *seq = new OptTaskSequence();
    *seq << "split:ho" << "paramsel:hoprimal" << "optimizer:rlsprimalrecinit";

    GurlsOptionsList * process = new GurlsOptionsList("processes", false);

    OptProcess* process1 = new OptProcess();
    *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
    process->addOpt("one", process1);

    this->opt->addOpt("seq", seq);
    this->opt->addOpt("processes", process);

    GURLS G;
    G.run(X, y, *(this->opt), std::string("one"));

    //const gMat2D<T>& times = this->opt->getOptValue<OptMatrix<gMat2D<T> > >("time.one");
    //std::cout << "train: " << times.getData()[2]*1000 << std::endl;
}

template <typename T>
void RecursiveRLSWrapper<T>::update(const gVec<T> &X, const gVec<T> &y)
{
    if(!this->trainedModel())
        throw gException("Error, Train Model First");

    RLSPrimalRecUpdate<T> optimizer;

    const unsigned long d = X.getSize();
    const unsigned long t = y.getSize();

    gMat2D<T>X_mat(1, d);
    copy(X_mat.getData(), X.getData(), d);
    gMat2D<T>y_mat(1, t);
    copy(y_mat.getData(), y.getData(), t);

    GurlsOptionsList* ret = optimizer.execute(X_mat, y_mat, *(this->opt));
    this->opt->removeOpt("optimizer");
    this->opt->addOpt("optimizer", ret);
}

template <typename T>
gMat2D<T>* RecursiveRLSWrapper<T>::eval(const gMat2D<T> &X)
{
    if(!this->trainedModel())
        throw gException("Error, Train Model First");

    gurls::PredPrimal<T> predTask;
    gMat2D<T> y;

    OptMatrix<gMat2D<T> >* result = predTask.execute(X, y, *(this->opt));

    gMat2D<T>& pred_mat = result->getValue();
    result->detachValue();
    delete result;

    return &pred_mat;
}

template <typename T>
RecursiveRLSRetrainWrapper<T>::RecursiveRLSRetrainWrapper(const std::string &name): RecursiveRLSWrapper<T>(name) {}

boost::posix_time::ptime begin, end;
boost::posix_time::time_duration diff;

template <typename T>
void RecursiveRLSRetrainWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    if(this->opt != NULL)
        delete this->opt;

    this->opt = new GurlsOptionsList(this->name, true);

    this->opt->removeOpt("nholdouts");
    this->opt->addOpt("nholdouts", new OptNumber(1));

    SplitHo<T> splitTask;
    GurlsOptionsList* split = splitTask.execute(X, y, *(this->opt));
    this->opt->addOpt("split", split);


    const gMat2D<unsigned long>& split_indices = split->getOptValue<OptMatrix<gMat2D<unsigned long> > >("indices");
    const gMat2D<unsigned long>& split_lasts = split->getOptValue<OptMatrix<gMat2D<unsigned long> > >("lasts");

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = y.cols();

    const unsigned long last = split_lasts.getData()[0];
    const unsigned long nva = n-last;

    unsigned long* va = new unsigned long[nva];
    copy(va, split_indices.getData()+last, nva);

    gMat2D<T>* Xva = new gMat2D<T>(nva, d);
    gMat2D<T>* yva = new gMat2D<T>(nva, t);

    subMatrixFromRows(X.getData(), n, d, va, nva, Xva->getData());
    subMatrixFromRows(y.getData(), n, t, va, nva, yva->getData());

    gMat2D<T>* XtX = new gMat2D<T>(d, d);
    gMat2D<T>* Xty = new gMat2D<T>(d, t);

    begin = boost::posix_time::microsec_clock::local_time();

    dot(X.getData(), X.getData(), XtX->getData(), n, d, n, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
    dot(X.getData(), y.getData(), Xty->getData(), n, d, n, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

    end = boost::posix_time::microsec_clock::local_time();
    diff = end-begin;
    std::cout << "mult: " << diff.total_milliseconds() << std::endl;

    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
    kernel->addOpt("XtX", new OptMatrix<gMat2D<T> >(*XtX));
    kernel->addOpt("Xty", new OptMatrix<gMat2D<T> >(*Xty));

    kernel->addOpt("Xva", new OptMatrix<gMat2D<T> >(*Xva));
    kernel->addOpt("yva", new OptMatrix<gMat2D<T> >(*yva));

    nTot = n;
    this->opt->addOpt("kernel", kernel);

    ParamSelHoPrimal<T> paramselTask;
    this->opt->addOpt("paramsel", paramselTask.execute(X, y, *(this->opt)));


    RLSPrimalRecInit<T> optimizerTask;
    begin = boost::posix_time::microsec_clock::local_time();
    this->opt->addOpt("optimizer", optimizerTask.execute(X, y, *(this->opt)));
    end = boost::posix_time::microsec_clock::local_time();
    diff = end-begin;

    std::cout << "train: " << diff.total_milliseconds() << std::endl;
}

template <typename T>
void RecursiveRLSRetrainWrapper<T>::update(const gVec<T> &X, const gVec<T> &y)
{
    RecursiveRLSWrapper<T>::update(X, y);

    const unsigned long d = X.getSize();
    const unsigned long t = y.getSize();

    ++nTot;

    gMat2D<T>* xtx = new gMat2D<T>(d,d);
    gMat2D<T>* xty = new gMat2D<T>(d,t);

    dot(X.getData(), X.getData(), xtx->getData(), 1, d, 1, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
    dot(X.getData(), y.getData(), xty->getData(), 1, d, 1, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

    GurlsOptionsList* kernel = this->opt->template getOptAs<GurlsOptionsList>("kernel");

    const gMat2D<T>& XtX = kernel->getOptValue<OptMatrix<gMat2D<T> > >("XtX");
    const gMat2D<T>& Xty = kernel->getOptValue<OptMatrix<gMat2D<T> > >("Xty");

    axpy(d*d, (T)1.0, XtX.getData(), 1, xtx->getData(), 1);
    axpy(d*t, (T)1.0, Xty.getData(), 1, xty->getData(), 1);

    kernel->removeOpt("XtX");
    kernel->addOpt("XtX", new OptMatrix<gMat2D<T> >(*xtx));

    kernel->removeOpt("Xty");
    kernel->addOpt("Xty", new OptMatrix<gMat2D<T> >(*xty));


    unsigned long proportion = static_cast<unsigned long>(gurls::round(1.0/this->opt->getOptAsNumber("hoproportion")));

    if(nTot % proportion == 0)
    {
        const gMat2D<T>& Xva = kernel->getOptValue<OptMatrix<gMat2D<T> > >("Xva");
        const gMat2D<T>& yva = kernel->getOptValue<OptMatrix<gMat2D<T> > >("yva");

        const unsigned long nva = Xva.rows();
        const unsigned long nva_new = nva+1;


        gMat2D<T>* Xva_new = new gMat2D<T>(nva_new, d);

        const T* old_it = Xva.getData();
        T* new_it = Xva_new->getData();
        for(const T* end = new_it+(nva_new*d); new_it< end; old_it+=nva, new_it +=nva_new)
            copy(new_it, old_it, nva);

        copy(Xva_new->getData()+nva, X.getData(), d, nva_new, 1);

        kernel->removeOpt("Xva");
        kernel->addOpt("Xva", new OptMatrix<gMat2D<T> >(*Xva_new));


        gMat2D<T>* yva_new = new gMat2D<T>(nva_new, t);

        old_it = yva.getData();
        new_it = yva_new->getData();
        for(const T* end = new_it+(nva_new*t); new_it< end; old_it+=nva, new_it +=nva_new)
            copy(new_it, old_it, nva);

        copy(yva_new->getData()+nva, y.getData(), t, nva_new, 1);

        kernel->removeOpt("yva");
        kernel->addOpt("yva", new OptMatrix<gMat2D<T> >(*yva_new));

    }
}

template <typename T>
void RecursiveRLSRetrainWrapper<T>::retrain()
{
    GurlsOptionsList* kernel = this->opt->template getOptAs<GurlsOptionsList>("kernel");
    const gMat2D<T> &Xva = kernel->getOptValue<OptMatrix<gMat2D<T> > >("Xva");
    const gMat2D<T> &yva = kernel->getOptValue<OptMatrix<gMat2D<T> > >("yva");

    const unsigned long nva = Xva.rows();

    this->opt->removeOpt("paramsel");
    this->opt->removeOpt("optimizer");

    GurlsOptionsList* split = this->opt->template getOptAs<GurlsOptionsList>("split");
    split->removeOpt("indices");
    split->removeOpt("lasts");


    gMat2D<unsigned long>* indices = new gMat2D<unsigned long>(1, nTot);
    gMat2D<unsigned long>* lasts = new gMat2D<unsigned long>(1, 1);

    set(indices->getData(), 0ul, nTot-nva);
    unsigned long * it = indices->getData() + (nTot-nva);
    for(unsigned long i=0; i<nva; ++i, ++it)
        *it = i;

    lasts->getData()[0] = (nTot-nva);

    split->addOpt("indices", new OptMatrix<gMat2D<unsigned long> >(*indices));
    split->addOpt("lasts", new OptMatrix<gMat2D<unsigned long> >(*lasts));


    ParamSelHoPrimal<T> paramselTask;
    this->opt->addOpt("paramsel", paramselTask.execute(Xva, yva, *(this->opt)));


    RLSPrimalRecInit<T> optimizerTask;
    gMat2D<T> emptyMat;
    this->opt->addOpt("nTot", new OptNumber(nTot));
    this->opt->addOpt("optimizer", optimizerTask.execute(emptyMat, emptyMat, *(this->opt)));
    this->opt->removeOpt("nTot");

}

}

