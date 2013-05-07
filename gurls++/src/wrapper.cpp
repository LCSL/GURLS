#include "wrapper.h"

#include "gurls.h"

using namespace gurls;


GurlsWrapper::GurlsWrapper(const std::string& name):opt(NULL), name(name) {}

GurlsWrapper::~GurlsWrapper()
{
    if(opt != NULL)
        delete opt;
}

GurlsWrapper::T GurlsWrapper::eval(const gVec<T> &X, unsigned long *index)
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

const GurlsOptionsList &GurlsWrapper::getOpt() const
{
    return *opt;
}

bool GurlsWrapper::trainedModel()
{
    if(opt == NULL)
        return false;

    return opt->hasOpt("optimizer");
}



RecursiveRLSWrapper::RecursiveRLSWrapper(const std::string& name):GurlsWrapper(name) {}

void RecursiveRLSWrapper::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    if(opt != NULL)
        delete opt;

    opt = new GurlsOptionsList(name, true);

    OptTaskSequence *seq = new OptTaskSequence();
    *seq << "split:ho" << "paramsel:hoprimal" << "optimizer:rlsprimalrecinit";

    GurlsOptionsList * process = new GurlsOptionsList("processes", false);

    OptProcess* process1 = new OptProcess();
    *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
    process->addOpt("one", process1);

    opt->addOpt("seq", seq);
    opt->addOpt("processes", process);

    GURLS G;
    G.run(X, y, *opt, std::string("one"));
}

void RecursiveRLSWrapper::update(const gVec<T> &X, const gVec<T> &y)
{
    if(!trainedModel())
        throw gException("Error, Train Model First");

    RLSPrimalRecUpdate<T> optimizer;

    const unsigned long d = X.getSize();
    const unsigned long t = y.getSize();

    gMat2D<T>X_mat(1, d);
    copy(X_mat.getData(), X.getData(), d);
    gMat2D<T>y_mat(1, t);
    copy(y_mat.getData(), y.getData(), t);

    GurlsOptionsList* ret = optimizer.execute(X_mat, y_mat, *opt);
    opt->removeOpt("optimizer");
    opt->addOpt("optimizer", ret);
}

gMat2D<GurlsWrapper::T>* RecursiveRLSWrapper::eval(const gMat2D<T> &X)
{
    if(!trainedModel())
        throw gException("Error, Train Model First");

    gurls::PredPrimal<T> predTask;
    gMat2D<T> y;

    OptMatrix<gMat2D<T> >* result = predTask.execute(X, y, *opt);

    gMat2D<T>& pred_mat = result->getValue();
    result->detachValue();
    delete result;

    return &pred_mat;
}


RecursiveRLSRetrainWrapper::RecursiveRLSRetrainWrapper(const std::string &name): RecursiveRLSWrapper(name) {}

void RecursiveRLSRetrainWrapper::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    if(opt != NULL)
        delete opt;

    opt = new GurlsOptionsList(name, true);

    opt->removeOpt("nholdouts");
    opt->addOpt("nholdouts", new OptNumber(1));

    SplitHo<T> splitTask;
    GurlsOptionsList* split = splitTask.execute(X, y, *opt);
    opt->addOpt("split", split);


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

    dot(X.getData(), X.getData(), XtX->getData(), n, d, n, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
    dot(X.getData(), y.getData(), Xty->getData(), n, d, n, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
    kernel->addOpt("XtX", new OptMatrix<gMat2D<T> >(*XtX));
    kernel->addOpt("Xty", new OptMatrix<gMat2D<T> >(*Xty));

    kernel->addOpt("Xva", new OptMatrix<gMat2D<T> >(*Xva));
    kernel->addOpt("yva", new OptMatrix<gMat2D<T> >(*yva));

    nTot = n;
    opt->addOpt("kernel", kernel);

    ParamSelHoPrimal<T> paramselTask;
    opt->addOpt("paramsel", paramselTask.execute(X, y, *opt));

    RLSPrimalRecInit<T> optimizerTask;
    opt->addOpt("optimizer", optimizerTask.execute(X, y, *opt));
}

void RecursiveRLSRetrainWrapper::update(const gVec<GurlsWrapper::T> &X, const gVec<GurlsWrapper::T> &y)
{
    RecursiveRLSWrapper::update(X, y);

    const unsigned long d = X.getSize();
    const unsigned long t = y.getSize();

    ++nTot;

    gMat2D<T>* xtx = new gMat2D<T>(d,d);
    gMat2D<T>* xty = new gMat2D<T>(d,t);

    dot(X.getData(), X.getData(), xtx->getData(), 1, d, 1, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
    dot(X.getData(), y.getData(), xty->getData(), 1, d, 1, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

    GurlsOptionsList* kernel = opt->getOptAs<GurlsOptionsList>("kernel");

    const gMat2D<T>& XtX = kernel->getOptValue<OptMatrix<gMat2D<T> > >("XtX");
    const gMat2D<T>& Xty = kernel->getOptValue<OptMatrix<gMat2D<T> > >("Xty");

    axpy(d*d, (T)1.0, XtX.getData(), 1, xtx->getData(), 1);
    axpy(d*t, (T)1.0, Xty.getData(), 1, xty->getData(), 1);

    kernel->removeOpt("XtX");
    kernel->addOpt("XtX", new OptMatrix<gMat2D<T> >(*xtx));

    kernel->removeOpt("Xty");
    kernel->addOpt("Xty", new OptMatrix<gMat2D<T> >(*xty));


    unsigned long proportion = static_cast<unsigned long>(gurls::round(1.0/opt->getOptAsNumber("hoproportion")));

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

void RecursiveRLSRetrainWrapper::retrain()
{
    GurlsOptionsList* kernel = opt->getOptAs<GurlsOptionsList>("kernel");
    const gMat2D<T> &Xva = kernel->getOptValue<OptMatrix<gMat2D<T> > >("Xva");
    const gMat2D<T> &yva = kernel->getOptValue<OptMatrix<gMat2D<T> > >("yva");

    const unsigned long nva = Xva.rows();

    opt->removeOpt("paramsel");
    opt->removeOpt("optimizer");

    GurlsOptionsList* split = opt->getOptAs<GurlsOptionsList>("split");
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
    opt->addOpt("paramsel", paramselTask.execute(Xva, yva, *opt));


    RLSPrimalRecInit<T> optimizerTask;
    gMat2D<T> emptyMat;
    opt->addOpt("nTot", new OptNumber(nTot));
    opt->addOpt("optimizer", optimizerTask.execute(emptyMat, emptyMat, *opt));
    opt->removeOpt("nTot");

}


