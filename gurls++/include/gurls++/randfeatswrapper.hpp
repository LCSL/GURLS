#include "gurls++/randfeatswrapper.h"

#include "gurls++/splitho.h"
#include "gurls++/hoprimal.h"
#include "gurls++/rlsprimal.h"
#include "gurls++/primal.h"
#include "gurls++/options.h"
#include "gurls++/utils.h"


namespace gurls
{

template <typename T>
RandomFeaturesWrapper<T>::RandomFeaturesWrapper(const std::string &name): RLSWrapper<T>(name), W(NULL) { }

template <typename T>
RandomFeaturesWrapper<T>::~RandomFeaturesWrapper()
{
    if(W != NULL)
        delete W;
}

template <typename T>
void RandomFeaturesWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
//    D = opt.n_randfeats;
    const unsigned long D = static_cast<unsigned long>(this->opt->getOptAsNumber("randfeats.D"));

//    [n,d] = size(Xtr);
    const unsigned long d = X.cols();


//    W = sqrt(2)*randn(d,D);

    if(W != NULL)
        delete W;

    W = rp_projections<T>(d, D);


//    V = X*W;
//    Xtr = [cos(V) sin(V)];
    gMat2D<T> *Xtr = rp_apply_real(X, *W);

    RLSWrapper<T>::train(*Xtr, y);

    delete Xtr;
}

template<typename T>
gMat2D<T>* RandomFeaturesWrapper<T>::eval(const gMat2D<T> &X)
{
    if(W == NULL)
        throw gException("Error, Train Model First");


//    V = X*W;
//    Xte = [cos(V) sin(V)];
    gMat2D<T> *Xte = rp_apply_real(X, *W);

    gMat2D<T>* pred = RLSWrapper<T>::eval(*Xte);

    delete Xte;

    return pred;
}

template<typename T>
void RandomFeaturesWrapper<T>::setNRandFeats(unsigned long value)
{
    this->opt->template getOptValue<OptNumber>("randfeats.D") = value;
}

}
