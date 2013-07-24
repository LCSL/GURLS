#include "gurls++/randfeatswrapper.h"

#include "gurls++/splitho.h"
#include "gurls++/hoprimal.h"
#include "gurls++/rlsprimal.h"
#include "gurls++/primal.h"
#include "gurls++/options.h"

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>


namespace gurls
{

template <typename T>
RandomFeaturesWrapper<T>::RandomFeaturesWrapper(const std::string &name): RLSWrapper<T>(name), W(NULL)
{
    this->opt->addOpt("n_randfeats", new OptNumber(100));
}

template <typename T>
void RandomFeaturesWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{

//    D = opt.n_randfeats;
    const unsigned long D = static_cast<unsigned long>(this->opt->getOptAsNumber("n_randfeats"));

//    [n,d] = size(Xtr);
    const unsigned long n = X.rows();
    const unsigned long d = X.cols();


//    W = sqrt(2)*randn(d,D);

    boost::random::mt19937 gen;
    boost::random::normal_distribution<> g(0.0, sqrt(2.0));

    if(W != NULL)
        delete W;

    W = new gMat2D<T>(d, D);
    for(T* W_it = W->getData(), *const W_end = W_it + (d*D); W_it != W_end; ++W_it)
        *W_it = g(gen);

//    V = X*W;
    T* V = new T[n*D];
    dot(X.getData(), W->getData(), V, n, d, d, D, n, D, CblasNoTrans, CblasNoTrans, CblasColMajor);


//    Xtr = [cos(V) sin(V)];
    gMat2D<T> Xtr (n, 2*D);
    T* Xtr_it = Xtr.getData();
    for(T *Xtr_end = Xtr_it + (n*D), *V_it = V; Xtr_it != Xtr_end; ++Xtr_it, ++V_it)
        *Xtr_it = cos(*V_it);
    for(T *Xtr_end = Xtr_it + (n*D), *V_it = V; Xtr_it != Xtr_end; ++Xtr_it, ++V_it)
        *Xtr_it = sin(*V_it);

    delete [] V;

    RLSWrapper<T>::train(Xtr, y);
}

template<typename T>
gMat2D<T>* RandomFeaturesWrapper<T>::eval(const gMat2D<T> &X)
{
    if(W == NULL)
        throw gException("Error, Train Model First");

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();

    const unsigned long D = static_cast<unsigned long>(this->opt->getOptAsNumber("n_randfeats"));

//    V = X*W;
    T* V = new T[n*D];
    dot(X.getData(), W->getData(), V, n, d, d, D, n, D, CblasNoTrans, CblasNoTrans, CblasColMajor);

//    Xte = [cos(V) sin(V)];
    gMat2D<T> Xte (n, 2*D);
    T* Xte_it = Xte.getData();
    for(T *Xte_end = Xte_it + (n*D), *V_it = V; Xte_it != Xte_end; ++Xte_it, ++V_it)
        *Xte_it = cos(*V_it);
    for(T *Xte_end = Xte_it + (n*D), *V_it = V; Xte_it != Xte_end; ++Xte_it, ++V_it)
        *Xte_it = sin(*V_it);

    delete [] V;

    return RLSWrapper<T>::eval(Xte);
}

template<typename T>
void RandomFeaturesWrapper<T>::setNRandFeats(unsigned long value)
{
    this->opt->template getOptValue<OptNumber>("n_randfeats") = value;
}

}
