
#include "gurls++/nystromwrapper.h"

#include "gurls++/optmatrix.h"
#include "gurls++/predkerneltraintest.h"
#include "gurls++/primal.h"
#include "gurls++/utils.h"


#include "gurls++/macroavg.h"
#include "gurls++/precisionrecall.h"
#include "gurls++/rmse.h"

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <vector>

namespace gurls
{

template <typename T>
NystromWrapper<T>::NystromWrapper(const std::string& name):KernelWrapper<T>(name) {}

template <typename T>
void NystromWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    GurlsOptionsList *opt = this->opt;
    const bool regression = (this->probType == GurlsWrapper<T>::REGRESSION);

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
    const unsigned long nparams = static_cast<unsigned long>(opt->getOptAsNumber("nlambda"));


    const gMat2D<T> *Xtr = NULL;
    const gMat2D<T> *ytr = NULL;
    const gMat2D<T> *Xva = NULL;
    const gMat2D<T> *yva = NULL;


    bool split = (nparams > 1);

    unsigned long ntr, nva;

    if(split)
    {
        ntr = static_cast<unsigned long>(n*(1.0-opt->getOptAsNumber("hoproportion")));
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

    std::vector<unsigned long> guesses;
    unsigned long step = static_cast<unsigned long>(floor(n_nystrom/nparams));
    for(unsigned long i = guesses_end-1, count = 0; count < nparams; i-=step, ++count)
        guesses.push_back(i);


//    indices = randperm(ntr);
    unsigned long indices_length = 0;
    unsigned long *indices = getIndices(y, ntr, t, guesses_end, indices_length);



    gMat2D<unsigned long> *indices_mat = new gMat2D<unsigned long>(indices_length, 1);
    copy(indices_mat->getData(), indices, indices_length);
    paramsel->addOpt("indices", new OptMatrix<gMat2D<unsigned long> >(*indices_mat));


//    vout.perf = zeros(length(guesses),1);
    gMat2D<T> *perf_mat = new gMat2D<T>(guesses.size(), 1);
    paramsel->addOpt("acc", new OptMatrix<gMat2D<T> >(*perf_mat));
    T* perf_it = perf_mat->getData();
    set(perf_it, (T)-1, guesses.size());

    gMat2D<unsigned long> *times_mat = new gMat2D<unsigned long>(guesses.size(), 1);
    paramsel->addOpt("times", new OptMatrix<gMat2D<unsigned long> >(*times_mat));
    unsigned long *times_it = times_mat->getData();
    set(times_it, 0ul, guesses.size());



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

    gMat2D<unsigned long> *guesses_mat = new gMat2D<unsigned long>(guesses.size(), 1);
    unsigned long *guesses_mat_it = guesses_mat->getData();
    set(guesses_mat_it, 0ul, guesses.size());

    boost::posix_time::ptime begin, end;
    boost::posix_time::time_duration diff;

    T prevPerf = (regression)? -std::numeric_limits<T>::max(): 0;

    begin = boost::posix_time::microsec_clock::local_time();
    for(std::vector<unsigned long>::iterator it = guesses.begin(); it != guesses.end(); ++it)
    {
        const unsigned long i = *it;

        *guesses_mat_it = i+1;

//        Xnew = Xtr(indices(i_init:i_end),:);
        const unsigned long nindices = i - i_init + 1;
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

            if(K_mat != NULL)
                delete K_mat;

            K_mat = &predkernel_K;
            K = K_mat->getData();
        }


//        KtK(i_init:i_end,i_init:i_end) = kernel.K'*kernel.K;
        T* KtK_sub = new T[nindices*nindices];

        dot(predkernel_K.getData(), predkernel_K.getData(), KtK_sub, ntr, nindices, ntr, nindices, nindices, nindices, CblasTrans, CblasNoTrans, CblasColMajor);
        for(unsigned long i=0; i< nindices; ++i)
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

            if(regression)
                *perf_it =  static_cast<T>(perf->getOptValue<OptNumber>("forho"));
            else
            {
        //        vout.perf(guesses_c) = mean(perf.forho);
                const gMat2D<T> &acc = perf->getOptValue<OptMatrix<gMat2D<T> > >("acc");

                *perf_it = sumv(acc.getData(), acc.getSize())/acc.getSize();
            }


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

        *times_it = static_cast<unsigned long>(diff.total_milliseconds());
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
    }

    delete K_mat;

    delete [] indices;
    delete [] KtK;

    delete perfTask;
    delete opt_tmp;


    unsigned long nsub;

    if(guesses_mat_it - guesses_mat->getData() < guesses_mat->getSize()-1)
    {
        nsub = *guesses_mat_it;
        gMat2D<T> *Xsub_mat_sub = new gMat2D<T>(nsub, d);

        const unsigned long size = nsub*d;

        for(T* it = Xsub_mat_sub->getData(), *end = it+(size), *X_it = Xsub_mat->getData(); it != end; it+=nsub, X_it+=guesses_end)
            copy(it, X_it, nsub);

        delete Xsub_mat;
        Xsub_mat = Xsub_mat_sub;
    }
    else
        nsub = guesses_end;


//    vout.X = Xsub(1:i_end,:);
    paramsel->addOpt("X", new OptMatrix<gMat2D<T> >(*Xsub_mat));

//    vout.C = alpha;
    paramsel->addOpt("C", new OptMatrix<gMat2D<T> >(*alpha_mat));

//    [dummy,vout.n_opt] = max(vout.perf);
    paramsel->addOpt("n_opt", new OptNumber(*std::max_element(perf_mat->getData(), perf_mat->getData()+perf_mat->getSize())));

//    vout.guesses = guesses;
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<unsigned long> >(*guesses_mat));

}

template <typename T>
void NystromWrapper<T>::train_largescale(const gMat2D<T> &X, const gMat2D<T> &y)
{
    GurlsOptionsList *opt = this->opt;

    GurlsOptionsList* paramsel = opt->getOptAs<GurlsOptionsList>("paramsel");
    paramsel->removeOpt("X");
    paramsel->removeOpt("C");
    paramsel->removeOpt("n_opt");
    paramsel->removeOpt("guesses");
    paramsel->removeOpt("acc");
    paramsel->removeOpt("times");
    paramsel->removeOpt("indices");


    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = y.cols();

    const gMat2D<T> empty;

    const double sigma = opt->getOptAsNumber("paramsel.sigma");

    const double n_nystrom = opt->getOptAsNumber("n_nystrom");
//    B = opt.nystrom_blocksize;
    const unsigned long B = static_cast<unsigned long>(opt->getOptAsNumber("nlambda"));

    const gMat2D<T> *Xtr = &X;
    const gMat2D<T> *ytr = &y;

    unsigned long ntr = n;


    unsigned long guesses_end = static_cast<unsigned long>(n_nystrom);

//    guesses = B:B:opt.rank_max;
    std::vector<unsigned long> guesses;
    for(unsigned long i = B-1; i < guesses_end; i+=B)
        guesses.push_back(i);

    guesses_end = *(guesses.rbegin())+1;

//    indices = randperm(ntr);
    unsigned long indices_length = 0;
    unsigned long *indices = getIndices(y, ntr, t, guesses_end, indices_length);



    gMat2D<unsigned long> *indices_mat = new gMat2D<unsigned long>(indices_length, 1);
    copy(indices_mat->getData(), indices, indices_length);
    paramsel->addOpt("indices", new OptMatrix<gMat2D<unsigned long> >(*indices_mat));


    gMat2D<unsigned long> *times_mat = new gMat2D<unsigned long>(guesses.size(), 1);
    paramsel->addOpt("times", new OptMatrix<gMat2D<unsigned long> >(*times_mat));
    unsigned long *times_it = times_mat->getData();
    set(times_it, 0ul, guesses.size());



//    Xsub = zeros(guesses(end),d);
    gMat2D<T> *Xsub_mat = new gMat2D<T>(guesses_end, d);
    T *Xsub = Xsub_mat->getData();
    set(Xsub, (T)0.0, guesses_end*d);

//    KtK = zeros(guesses(end),guesses(end));
    T *KtK = new T[guesses_end*guesses_end];
    set(KtK, (T)0.0, guesses_end*guesses_end);

    //    Kty = zeros(guesses(end),T);
    T *Kty = new T[guesses_end*t];
    set(Kty, (T)0.0, guesses_end*t);

//    i_init = 1;
    unsigned long i_init = 0;


//    opt_tmp.paramsel.sigma = opt.paramsel.sigma;
//    opt_tmp.kernel.type = opt.kernel.type;

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

    gMat2D<unsigned long> *guesses_mat = new gMat2D<unsigned long>(guesses.size(), 1);
    unsigned long *guesses_mat_it = guesses_mat->getData();
    set(guesses_mat_it, 0ul, guesses.size());

    boost::posix_time::ptime begin, end;
    boost::posix_time::time_duration diff;

    begin = boost::posix_time::microsec_clock::local_time();
    for(std::vector<unsigned long>::iterator it = guesses.begin(); it != guesses.end(); ++it)
    {
        const unsigned long i = *it;
        *guesses_mat_it = i+1;

//        Xnew = Xtr(indices(i_init:i_end),:);
        const unsigned long nindices = i - i_init + 1;
        gMat2D<T> *Xnew = new gMat2D<T>(nindices,d);
        subMatrixFromRows(Xtr->getData(), ntr, d, indices+i_init, nindices, Xnew->getData());

//        opt_tmp.rls.X = Xnew;
        tmp_optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*Xnew, false));

//        kernel_col = predkernel_traintest(Xtr,[],opt_tmp);
        GurlsOptionsList *kernel = predKernelTask.execute(*Xtr, empty, *opt_tmp);

        gMat2D<T> &predkernel_K = kernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

//        Kty(i_init:i_end,:) = kernel_col.K'*y;
        T* Kty_sub = new T[nindices*t];
        dot(predkernel_K.getData(), ytr->getData(), Kty_sub, ntr, nindices, ntr, t, nindices, t, CblasTrans, CblasNoTrans, CblasColMajor);
        for(unsigned long i=0; i< t; ++i)
            copy(Kty+i_init+(i*guesses_end), Kty_sub+(i*nindices), nindices);

        delete [] Kty_sub;

//        KtK(i_init:i_end,i_init:i_end) = kernel_col.K'*kernel_col.K;
        T* KtK_sub = new T[nindices*nindices];

        dot(predkernel_K.getData(), predkernel_K.getData(), KtK_sub, ntr, nindices, ntr, nindices, nindices, nindices, CblasTrans, CblasNoTrans, CblasColMajor);
        for(unsigned long i=0; i< nindices; ++i)
            copy(KtK+i_init+((i_init+i)*guesses_end), KtK_sub+(i*nindices), nindices);

        delete[] KtK_sub;


//        if guesses_c>1
        if(it != guesses.begin())
        {
//            KtKcol = zeros((i_init-1),i_end-i_init+1);
            T *KtKcol = new T[i_init*nindices];
            set(KtKcol, (T)0.0, i_init*nindices);


            const T salpha = (T)(-1.0/pow(sigma, 2));
            const T one = (T) 1.0;

            T *kernel_old = new T[i_init];

            const T *Xtr_it = Xtr->getData();
            T *K_it = predkernel_K.getData();

//            for l=1:ntot
            for(unsigned long l = 0; l< ntr; ++l, ++Xtr_it, ++K_it)
            {
//                opt_tmp.rls.X = Xsub(1:(i_init-1),:);
//                kernel_old = predkernel_traintest(X(l,:),[],opt_tmp);
                distance_transposed_vm(Xtr_it, Xsub, d, guesses_end, kernel_old, i_init, ntr);
                scal(i_init, salpha, kernel_old, 1);
                exp(kernel_old, i_init);

//                KtKcol = KtKcol + kernel_old.K'*kernel_col.K(l,:);
                gemv(CblasNoTrans, i_init, 1, one, kernel_old, i_init, K_it, ntr, one, KtKcol, 1); // 1 i_init   1 nindices  i_init nindices
            }

            delete [] kernel_old;


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


//        alpha = pinv(KtK(1:i_end,1:i_end))*(Kty(1:i_end,:));

        const unsigned long ii = i+1;

        KtK_sub = new T[ ii*ii];
        for(T *K_it = KtK, *Ks_it = KtK_sub, *const K_end = K_it+(guesses_end*ii); K_it != K_end; K_it+=guesses_end, Ks_it+=ii)
            copy(Ks_it, K_it, ii);

        int r, c;
        T *pinv_K = pinv(KtK_sub, ii, ii, r, c);
        delete [] KtK_sub;

        Kty_sub = new T[ii*t];
        for(T *K_it = Kty, *Ks_it = Kty_sub, *const K_end = K_it+(guesses_end*t); K_it != K_end; K_it+=guesses_end, Ks_it+=ii)
            copy(Ks_it, K_it, ii);


        dot(pinv_K, Kty_sub, alpha, ii, ii, ii, t, ii, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

        delete [] Kty_sub;
        delete [] pinv_K;


        tmp_optimizer->removeOpt("X");

//        Xsub(i_init:i_end,:) = Xnew;
        for(T* Xsub_it = Xsub+i_init, *Xsub_end = Xsub_it+(guesses_end*d), *Xnew_it = Xnew->getData(); Xsub_it != Xsub_end; Xsub_it+=guesses_end, Xnew_it+=nindices)
            copy(Xsub_it, Xnew_it, nindices);

        delete Xnew;

//        i_init = i_end+1;
        i_init = i+1;

        end = boost::posix_time::microsec_clock::local_time();
        diff = end-begin;

        *times_it = static_cast<unsigned long>(diff.total_milliseconds());
        ++times_it;

        ++guesses_mat_it;
    }


    delete [] indices;
    delete [] KtK;
    delete [] Kty;

    delete opt_tmp;


//    vout.X = Xsub(1:i_end,:);
    paramsel->addOpt("X", new OptMatrix<gMat2D<T> >(*Xsub_mat));

//    vout.C = alpha;
    paramsel->addOpt("C", new OptMatrix<gMat2D<T> >(*alpha_mat));

//    vout.guesses = guesses;
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<unsigned long> >(*guesses_mat));

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
gMat2D<T>* NystromWrapper<T>::eval_largescale(const gMat2D<T> &X)
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
void NystromWrapper<T>::setParam(double value)
{
    this->setNparams(1);

    if(this->opt->hasOpt("n_nystrom"))
        this->opt->template getOptValue<OptNumber>("n_nystrom") = value;
    else
        this->opt->addOpt("n_nystrom", new OptNumber(value));
}

//template <typename T>
//void NystromWrapper<T>::rescale(gMat2D<T> &y)
//{
//    const unsigned long n = y.rows();
//    const unsigned long t = y.cols();

//    T* yt = new T[y.getSize()];

//    transpose(y.getData(), n, t, yt);


//    //unsigned long counts[2] = {0, 0};
//    unsigned long *counts = new unsigned long [t];
//    set(counts, 0ul, t);

//    unsigned long *maxIndices = new unsigned long [n];
//    unsigned long *m_it = maxIndices;


//    T *it = yt;
//    for(unsigned long i=0; i<n; ++i, it+=t)
//    {
//        unsigned long mindex = std::max_element(it, it+t)-it;

//        *m_it++ = mindex;
//        ++counts[mindex];
//    }

//    delete [] yt;

//    unsigned long *coeffs = new unsigned long [t];
//    set(coeffs, 0ul, t);
//    for(unsigned long i=0; i<t; ++i)
//    {
//        for(unsigned long j=0; j<t; ++j)
//            if(i!=j)
//                coeffs[i] += counts[j];
//    }

//    it = y.getData();
//    m_it = maxIndices;

//    for(unsigned long i=0; i<n; ++i, ++it, ++m_it)
////        scal(t, (T)(counts[!(*m_it)]), it, n);
//        scal(t, (T)(coeffs[*m_it]), it, n);


//    delete [] counts;
//    delete [] coeffs;


//    delete [] maxIndices;

//}

template <typename T>
unsigned long* NystromWrapper<T>::getIndices(const gMat2D<T>&y, const unsigned long n, const unsigned long t, const unsigned long n_nystrom, unsigned long &length)
{
    unsigned long *indices = new unsigned long[n];
    randperm(n, indices, true, 0);

    length = n;
    return indices;
}

}
