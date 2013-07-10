
//#include "gurls.h"
//#include "optlist.h"
#include "recrlswrapper.h"

namespace gurls
{
template <typename T>
RecursiveRLSWrapper<T>::RecursiveRLSWrapper(const std::string &name): GurlsWrapper<T>(name)
{
    this->opt->template getOptValue<OptNumber>("nholdouts") = 1.0;
}

template <typename T>
void RecursiveRLSWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    this->opt->removeOpt("split");
    this->opt->removeOpt("paramsel");
    this->opt->removeOpt("optimizer");
    this->opt->removeOpt("kernel");


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


    dot(X.getData(), X.getData(), XtX->getData(), n, d, n, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
    dot(X.getData(), y.getData(), Xty->getData(), n, d, n, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);


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
    this->opt->addOpt("optimizer", optimizerTask.execute(X, y, *(this->opt)));
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
void RecursiveRLSWrapper<T>::retrain()
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
