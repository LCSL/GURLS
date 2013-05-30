
#include "nystromwrapper.h"

#include "optmatrix.h"
#include "predkerneltraintest.h"
#include "primal.h"
#include "utils.h"
//#include "splitho.h"

#include "macroavg.h"
#include "precisionrecall.h"
#include "rmse.h"

#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace gurls
{

template <typename T>
NystromWrapper<T>::NystromWrapper(const std::string& name):GurlsWrapper<T>(name)
{
    this->opt = new GurlsOptionsList(name, true);

    GurlsOptionsList *paramsel = new GurlsOptionsList("paramsel");
    this->opt->addOpt("paramsel", paramsel);

    GurlsOptionsList *kernel = new GurlsOptionsList("kernel");
    this->opt->addOpt("kernel", kernel);


    this->opt->removeOpt("hoproportion");
    this->opt->addOpt("hoproportion", new OptNumber(0.2));

    this->opt->removeOpt("nholdouts");
    this->opt->addOpt("nholdouts", new OptNumber(1));
}

template <typename T>
void NystromWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    GurlsOptionsList *opt = this->opt;

    GurlsOptionsList* paramsel = opt->getOptAs<GurlsOptionsList>("paramsel");
    paramsel->removeOpt("X");
    paramsel->removeOpt("C");
    paramsel->removeOpt("n_opt");
    paramsel->removeOpt("guesses");
    paramsel->removeOpt("perf");
    paramsel->removeOpt("perf_train");


    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = y.cols();

    const gMat2D<T> empty;

    const double sigma = opt->getOptAsNumber("paramsel.sigma");
    Performance<T> *perfTask = Performance<T>::factory(opt->getOptAsString("hoperf"));


    const double n_nystrom = opt->getOptAsNumber("n_nystrom");
    const unsigned long nparams = static_cast<unsigned long>(opt->getOptAsNumber("nparams"));


    const gMat2D<T> *Xtr = NULL;
    const gMat2D<T> *ytr = NULL;
    const gMat2D<T> *Xva = NULL;
    const gMat2D<T> *yva = NULL;

    bool split = (nparams > 1);

    unsigned long ntr, nva;

    if(split)
    {
//        SplitHo<T> splitTask;

//        GurlsOptionsList* splitOpt = splitTask.execute(empty, y, *opt);
//        const gMat2D<unsigned long>& split_indices = splitOpt->getOptValue<OptMatrix<gMat2D<unsigned long> > >("indices");
//        const gMat2D<unsigned long>& split_lasts = splitOpt->getOptValue<OptMatrix<gMat2D<unsigned long> > >("lasts");

//        const unsigned long ntr = split_lasts.getData()[0];
//        const unsigned long nva = split_indices.getSize()-ntr;

//        const unsigned long *tr = split_indices.getData();
//        const unsigned long *va = split_indices.getData()+ntr;

        ntr = n*(1.0-opt->getOptAsNumber("hoproportion"));
        nva = n-ntr;

        unsigned long* trva = new unsigned long[n];
        randperm(n, trva, true, 0);

        unsigned long *tr = trva;
        unsigned long *va = trva+ntr;


        gMat2D<T> *Xtr_tmp = new gMat2D<T>(ntr, d);
        subMatrixFromRows(X.getData(), n, d, tr, ntr, Xtr_tmp->getData());
        Xtr = Xtr_tmp;

        gMat2D<T> *ytr_tmp = new gMat2D<T> (ntr, t);
        subMatrixFromRows(y.getData(), n, t, tr, ntr, ytr_tmp->getData());
        ytr = ytr_tmp;

        gMat2D<T> *Xva_tmp = new gMat2D<T> (nva, d);
        subMatrixFromRows(X.getData(), n, d, va, nva, Xva_tmp->getData());
        Xva = Xva_tmp;

        gMat2D<T> *yva_tmp = new gMat2D<T> (nva, t);
        subMatrixFromRows(y.getData(), n, t, va, nva, yva_tmp->getData());
        yva = yva_tmp;


    //    delete splitOpt;
        delete [] trva;
    }
    else
    {
        Xtr = &X;
        ytr = &y;

        ntr = n;
        nva = 0;
    }



    unsigned long guesses_end = static_cast<unsigned long>(n_nystrom);

//    guesses = unique(round((opt.n_nystrom).^((1:opt.nparams)./opt.nparams)));
    std::set<unsigned long> guesses;
//    for(unsigned long i=1; i<=nparams; ++i)
//        guesses.insert( static_cast<unsigned long>(gurls::round(std::pow(n_nystrom, ((T)i)/nparams)))-1);

    unsigned long aaa = static_cast<unsigned long>(floor(n_nystrom/nparams));
    for(unsigned long i = guesses_end-1, count = 0; count < nparams; i-=aaa, ++count)
        guesses.insert(i);


//    indices = randperm(ntr);
    unsigned long *indices = new unsigned long[ntr];
    randperm(ntr, indices, true, 0);



//    vout.perf = zeros(length(guesses),1);
    gMat2D<T> *perf_mat = new gMat2D<T>(guesses.size(), 1);
    paramsel->addOpt("acc", new OptMatrix<gMat2D<T> >(*perf_mat));
    T* perf_it = perf_mat->getData();
    set(perf_it, (T)-1, guesses.size());

    gMat2D<T> *times_mat = new gMat2D<T>(guesses.size(), 1);
    paramsel->addOpt("times", new OptMatrix<gMat2D<T> >(*times_mat));
    T *times_it = times_mat->getData();
    set(times_it, (T)0, guesses.size());



//    Xsub = zeros(guesses(end),d);
    gMat2D<T> *Xsub_mat = new gMat2D<T>(guesses_end, d);
    T *Xsub = Xsub_mat->getData();
    set(Xsub, (T)0.0, guesses_end*d);

//    K = zeros(ntr,guesses(end));
//    Ktrva = zeros(nva,guesses(end));
    gMat2D<T>* K_mat = NULL;
    T *K = NULL;
    T* Ktrva = NULL;
    if(split)
    {
        K_mat = new gMat2D<T>(ntr, guesses_end);
        K = K_mat->getData();
        set(K, (T)0.0, ntr*guesses_end);

        Ktrva = new T[nva*guesses_end];
        set(Ktrva, (T)0.0, nva*guesses_end);
    }

//    KtK = zeros(guesses(end),guesses(end));
    T* KtK = new T[guesses_end*guesses_end];
    set(KtK, (T)0.0, guesses_end*guesses_end);

//    i_init = 1;
    unsigned long i_init = 0;


//    opt_tmp.paramsel.sigma = opt.paramsel.sigma;
    GurlsOptionsList* opt_tmp = new GurlsOptionsList("tmp");
    GurlsOptionsList* tmp_kernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* tmp_paramsel = new GurlsOptionsList("paramsel");
    GurlsOptionsList* tmp_optimizer = new GurlsOptionsList("optimizer");

    opt_tmp->addOpt("kernel", tmp_kernel);
    opt_tmp->addOpt("paramsel", tmp_paramsel);
    opt_tmp->addOpt("optimizer", tmp_optimizer);

    tmp_kernel->addOpt("type", opt->getOptAsString("kernel.type"));
    tmp_paramsel->addOpt("sigma", new OptNumber(sigma));


    PredKernelTrainTest<T> predKernelTask;

    gMat2D<T> *alpha_mat = new gMat2D<T>(guesses_end, t);
    T *const alpha = alpha_mat->getData();

    gMat2D<T> *guesses_mat = new gMat2D<T>(guesses.size(), 1);
    T *guesses_mat_it = guesses_mat->getData();
    set(guesses_mat_it, (T)-1, guesses.size());

    boost::posix_time::ptime begin, end;
    boost::posix_time::time_duration diff;

    T prevPerf = 0;

    begin = boost::posix_time::microsec_clock::local_time();
    for(std::set<unsigned long>::iterator it = guesses.begin(); it != guesses.end(); ++it)
    {
        const unsigned long i = *it;

        *guesses_mat_it = i+1;

//        Xnew = Xtr(indices(i_init:i_end),:);
        const int nindices = i - i_init + 1;
        gMat2D<T> *Xnew = new gMat2D<T>(nindices,d);
        subMatrixFromRows(Xtr->getData(), ntr, d, indices+i_init, nindices, Xnew->getData());

//        opt_tmp.rls.X = Xnew;
        tmp_optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*Xnew, false));

//        kernel = predkernel_traintest(Xtr,[],opt_tmp);
        GurlsOptionsList *kernel = predKernelTask.execute(*Xtr, empty, *opt_tmp);

        gMat2D<T> &predkernel_K = kernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

        if(split)
        {
//            K(:,i_init:i_end) = kernel.K;
            copy(K+(ntr*i_init), predkernel_K.getData(), predkernel_K.getSize());
        }
        else
        {
            OptMatrix<gMat2D<T> > * kk = kernel->getOptAs<OptMatrix<gMat2D<T> > >("K");
            kk->detachValue();
            K_mat = &predkernel_K;
            K = K_mat->getData();
        }


//        KtK(i_init:i_end,i_init:i_end) = kernel.K'*kernel.K;
        T* KtK_sub = new T[nindices*nindices];

        dot(predkernel_K.getData(), predkernel_K.getData(), KtK_sub, ntr, nindices, ntr, nindices, nindices, nindices, CblasTrans, CblasNoTrans, CblasColMajor);
        for(int i=0; i< nindices; ++i)
            copy(KtK+i_init+((i_init+i)*guesses_end), KtK_sub+(i*nindices), nindices);

        delete[] KtK_sub;


//        if guesses_c>1
        if(it != guesses.begin())
        {
//            KtKcol = K(:,1:(i_init-1))'*kernel.K;
            T *KtKcol = new T[i_init*nindices];
            dot(K, predkernel_K.getData(), KtKcol, ntr, i_init, ntr, nindices, i_init, nindices, CblasTrans, CblasNoTrans, CblasColMajor);

//            KtK(1:(i_init-1),i_init:i_end) = KtKcol;
//            KtK(i_init:i_end,1:(i_init-1)) = KtKcol';
            for(unsigned long j=0; j<nindices; ++j)
            {
                copy(KtK+(guesses_end*(i_init+j)), KtKcol+(i_init*j), i_init);
                copy(KtK+i_init+j, KtKcol+(i_init*j), i_init, guesses_end, 1);
            }

            delete [] KtKcol;
        }

        delete kernel;

//        alpha = pinv(KtK(1:i_end,1:i_end))*(K(:,1:i_end)'*ytr);

        const unsigned long ii = i+1;

        KtK_sub = new T[ ii*ii];
        for(T *K_it = KtK, *Ks_it = KtK_sub, *const K_end = K_it+(guesses_end*ii); K_it != K_end; K_it+=guesses_end, Ks_it+=ii)
            copy(Ks_it, K_it, ii);

        int r, c;
        T *pinv_K = pinv(KtK_sub, ii, ii, r, c);
        delete [] KtK_sub;

        T *Kty = new T[ii*t];
        dot(K, ytr->getData(), Kty, ntr, ii, ntr, t, ii, t, CblasTrans, CblasNoTrans, CblasColMajor);


        dot(pinv_K, Kty, alpha, ii, ii, ii, t, ii, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

        delete [] Kty;
        delete [] pinv_K;


        if(split)
        {
    //        validation error
    //        kernel = predkernel_traintest(Xva,[],opt_tmp);
            kernel = predKernelTask.execute(*Xva, empty, *opt_tmp);

    //        Ktrva(:,i_init:i_end) = kernel.K;
            const gMat2D<T> &predkernel_Kva = kernel->getOptValue<OptMatrix<gMat2D<T> > >("K");
            copy(Ktrva+(nva*i_init), predkernel_Kva.getData(), predkernel_Kva.getSize());
            delete kernel;


    //        opt.pred = Ktrva(:,1:i_end)*alpha;
            gMat2D<T>* pred = new gMat2D<T>(nva ,t);
            dot(Ktrva, alpha, pred->getData(), nva, ii, ii, t, nva, t, CblasNoTrans, CblasNoTrans, CblasColMajor);
            opt->addOpt("pred", new OptMatrix<gMat2D<T> >(*pred));

    //        perf = opt.hoperf([],yva,opt);
            GurlsOptionsList *perf = perfTask->execute(empty, *yva, *opt);
            opt->removeOpt("pred");

    //        vout.perf(guesses_c) = mean(perf.forho);
            const gMat2D<T> &acc = perf->getOptValue<OptMatrix<gMat2D<T> > >("acc");

            *perf_it = sumv(acc.getData(), acc.getSize())/acc.getSize();

            delete perf;
        }

        tmp_optimizer->removeOpt("X");

//        Xsub(i_init:i_end,:) = Xnew;
        for(T* Xsub_it = Xsub+i_init, *Xsub_end = Xsub_it+(guesses_end*d), *Xnew_it = Xnew->getData(); Xsub_it != Xsub_end; Xsub_it+=guesses_end, Xnew_it+=nindices)
            copy(Xsub_it, Xnew_it, nindices);

        delete Xnew;

//        i_init = i_end+1;
        i_init = i+1;

        end = boost::posix_time::microsec_clock::local_time();
        diff = end-begin;

        *times_it = diff.total_milliseconds();
        ++times_it;

        if(le(*perf_it, prevPerf) && i > 10)
            break;
        else
        {
            prevPerf = *perf_it;

            ++perf_it;
            ++guesses_mat_it;
        }
    }

    if(split)
    {
        delete Xtr;
        delete ytr;
        delete Xva;
        delete yva;

        delete [] Ktrva;
        delete K_mat;
    }

    delete [] indices;
    delete [] KtK;

    delete perfTask;
    delete opt_tmp;

    if(guesses_mat_it - guesses_mat->getData() < guesses_mat->getSize()-1)
    {
        const unsigned long nsub = static_cast<unsigned long>(*guesses_mat_it);
        gMat2D<T> *Xsub_mat_sub = new gMat2D<T>(nsub, d);

        const unsigned long size = nsub*d;

        for(T* it = Xsub_mat_sub->getData(), *end = it+(size), *X_it = Xsub_mat->getData(); it != end; it+=nsub, X_it+=guesses_end)
            copy(it, X_it, nsub);

        delete Xsub_mat;
        Xsub_mat = Xsub_mat_sub;
    }


//    vout.X = Xsub(1:i_end,:);
    paramsel->addOpt("X", new OptMatrix<gMat2D<T> >(*Xsub_mat));

//    vout.C = alpha;
    paramsel->addOpt("C", new OptMatrix<gMat2D<T> >(*alpha_mat));

//    [dummy,vout.n_opt] = max(vout.perf);
    paramsel->addOpt("n_opt", new OptNumber(*std::max_element(perf_mat->getData(), perf_mat->getData()+perf_mat->getSize())));

//    vout.guesses = guesses;
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

}

template <typename T>
void NystromWrapper<T>::update(const gVec<T> &X, const gVec<T> &y)
{

}

template <typename T>
gMat2D<T>* NystromWrapper<T>::eval(const gMat2D<T> &X)
{
    GurlsOptionsList *opt = this->opt;
    const gMat2D<T> &alpha_mat = opt->getOptValue<OptMatrix<gMat2D<T> > >("paramsel.C");
    const gMat2D<T> &X_mat = opt->getOptValue<OptMatrix<gMat2D<T> > >("paramsel.X");
    const T *const alpha = alpha_mat.getData();

    const unsigned long n = X.rows();
    const unsigned long t = alpha_mat.cols();
    const unsigned long n_nystrom = X_mat.rows();

    PredKernelTrainTest<T> predKernelTask;

    gMat2D<T> empty;


    GurlsOptionsList* opt_tmp = new GurlsOptionsList("tmp");
    GurlsOptionsList* tmp_kernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* tmp_paramsel = new GurlsOptionsList("paramsel");
    GurlsOptionsList* tmp_optimizer = new GurlsOptionsList("optimizer");

    opt_tmp->addOpt("kernel", tmp_kernel);
    opt_tmp->addOpt("paramsel", tmp_paramsel);
    opt_tmp->addOpt("optimizer", tmp_optimizer);

    tmp_kernel->addOpt("type", opt->getOptAsString("kernel.type"));
    tmp_paramsel->addOpt("sigma", new OptNumber(opt->getOptAsNumber("paramsel.sigma")));
    tmp_optimizer->addOpt("X", opt->getOpt("paramsel.X"));

    gMat2D<T> *y = new gMat2D<T>(n, t);

    GurlsOptionsList *predKernel = predKernelTask.execute(X, empty, *opt_tmp);
    const gMat2D<T> &predKernel_K = predKernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

    dot(predKernel_K.getData(), alpha, y->getData(), n, n_nystrom, n_nystrom, t, n, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

    tmp_optimizer->removeOpt("X", false);
    delete opt_tmp;
    delete predKernel;

    return y;
}

template <typename T>
gMat2D<T>* NystromWrapper<T>::eval_ls(const gMat2D<T> &X)
{
    GurlsOptionsList *opt = this->opt;
    const gMat2D<T> &alpha_mat = opt->getOptValue<OptMatrix<gMat2D<T> > >("paramsel.alpha");
    const gMat2D<T> &X_mat = opt->getOptValue<OptMatrix<gMat2D<T> > >("paramsel.X");
    const T *const alpha = alpha_mat.getData();

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = alpha_mat.cols();
    const unsigned long n_nystrom = X_mat.rows();

    gMat2D<T> Xi_mat(1, d);
    T* yi = new T[t];
    T* const Xi = Xi_mat.getData();

    gMat2D<T> *y_mat = new gMat2D<T>(n, t);
    T* y_it = y_mat->getData();

    GurlsOptionsList* opt_tmp = new GurlsOptionsList("tmp");
    GurlsOptionsList* tmp_kernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* tmp_paramsel = new GurlsOptionsList("paramsel");
    GurlsOptionsList* tmp_optimizer = new GurlsOptionsList("optimizer");

    opt_tmp->addOpt("kernel", tmp_kernel);
    opt_tmp->addOpt("paramsel", tmp_paramsel);
    opt_tmp->addOpt("optimizer", tmp_optimizer);

    tmp_kernel->addOpt("type", opt->getOptAsString("kernel.type"));
    tmp_paramsel->addOpt("sigma", new OptNumber(opt->getOptAsNumber("sigma")));
    tmp_optimizer->addOpt("X", opt->getOpt("paramsel.X"));


    PredKernelTrainTest<T> predKernelTask;
    const gMat2D<T> empty;

    for(int i=0; i<n; ++i, ++y_it)
    {
        getRow(X.getData(), n, d, i, Xi);

        GurlsOptionsList *predKernel = predKernelTask.execute(Xi_mat, empty, *opt_tmp);
        const gMat2D<T> &predKernel_K = predKernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

        dot(predKernel_K.getData(), alpha, yi, 1, n_nystrom, n_nystrom, t, 1, t, CblasNoTrans, CblasNoTrans, CblasColMajor);
        copy(y_it, yi, t, n, 1);

        delete predKernel;
    }

    tmp_optimizer->removeOpt("X", false);
    delete opt_tmp;
    delete [] yi;

    return y_mat;
}

template <typename T>
void NystromWrapper<T>::setNNystrom(unsigned long n_nystrom)
{
    this->opt->removeOpt("n_nystrom");
    this->opt->addOpt("n_nystrom", new OptNumber(n_nystrom));
}

template <typename T>
void NystromWrapper<T>::setNParams(unsigned long nparams)
{
    this->opt->removeOpt("nparams");
    this->opt->addOpt("nparams", new OptNumber(nparams));
}

template <typename T>
void NystromWrapper<T>::setKernelType(const std::string &type)
{
    GurlsOptionsList* kernel = this->opt->template getOptAs<GurlsOptionsList>("kernel");
    kernel->removeOpt("type");
    kernel->addOpt("type", type);
}

template <typename T>
void NystromWrapper<T>::setSigma(double sigma)
{
    GurlsOptionsList* paramsel = this->opt->template getOptAs<GurlsOptionsList>("paramsel");
    paramsel->removeOpt("sigma");
    paramsel->addOpt("sigma", new OptNumber(sigma));
}

}
