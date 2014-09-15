
#include "gurls++/recrlswrapperchol.h"

namespace gurls
{
template <typename T>
RecursiveRLSCholUpdateWrapper<T>::RecursiveRLSCholUpdateWrapper(const std::string &name): GurlsWrapper<T>(name)
{
    this->opt->template getOptValue<OptNumber>("nholdouts") = 1.0;
}

template <typename T>
void RecursiveRLSCholUpdateWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
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


    RLSPrimalRecInitCholesky<T> optimizerTask;
    this->opt->addOpt("optimizer", optimizerTask.execute(X, y, *(this->opt)));
}

template <typename T>
void RecursiveRLSCholUpdateWrapper<T>::update(const gMat2D<T> &X, const gMat2D<T> &y)
{
    if(!this->trainedModel())
        throw gException("Error, Train Model First");

    RLSPrimalRecUpdateCholesky<T> optimizer;

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
gMat2D<T>* RecursiveRLSCholUpdateWrapper<T>::eval(const gMat2D<T> &X)
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
void RecursiveRLSCholUpdateWrapper<T>::retrain()
{}

}
