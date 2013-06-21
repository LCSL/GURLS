#include "icholwrapper.h"

#include "optmatrix.h"
#include "predkerneltraintest.h"
#include "primal.h"
#include "utils.h"
#include "splitho.h"

#include "macroavg.h"
#include "precisionrecall.h"
#include "rmse.h"

#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace gurls
{

template <typename T>
ICholWrapper<T>::ICholWrapper(const std::string& name):GurlsWrapper<T>(name)
{
    this->opt = new GurlsOptionsList(name, true);

    GurlsOptionsList *paramsel = new GurlsOptionsList("paramsel");
    this->opt->addOpt("paramsel", paramsel);

    GurlsOptionsList *split = new GurlsOptionsList("split");
    this->opt->addOpt("split", split);
}

template <typename T>
void ICholWrapper<T>::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    GurlsOptionsList*opt = this->opt;

    const unsigned long m = static_cast<unsigned long>(opt->getOptAsNumber("paramsel.rank_max"));
    const unsigned long n_rank = static_cast<unsigned long>(opt->getOptAsNumber("paramsel.n_rank"));

    const double sigma = opt->getOptAsNumber("paramsel.sigma");

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = y.cols();

    Performance<T> * perfTask = Performance<T>::factory(opt->getOptAsString("hoperf"));

    bool computePred = (opt->hasOpt("split.Xva") && opt->hasOpt("split.yva"));


    GurlsOptionsList* tmp_optimizer = new GurlsOptionsList("optimizer");
    OptMatrix<const gMat2D<T> > *X_opt = new OptMatrix<const gMat2D<T> >(X, false);
    tmp_optimizer->addOpt("X", X_opt);

    const gMat2D<T>* Xva = NULL;
    const gMat2D<T>* yva = NULL;
    gMat2D<T>* predKernel_K = NULL;
    const gMat2D<T> empty;

    if(computePred)
    {
        Xva = &(opt->getOptValue<OptMatrix<gMat2D<T> > >("split.Xva"));
        yva = &(opt->getOptValue<OptMatrix<gMat2D<T> > >("split.yva"));


        GurlsOptionsList* opt_tmp = new GurlsOptionsList("tmp");
        GurlsOptionsList* tmp_kernel = new GurlsOptionsList("kernel");
        GurlsOptionsList* tmp_paramsel = new GurlsOptionsList("paramsel");

        opt_tmp->addOpt("kernel", tmp_kernel);
        opt_tmp->addOpt("paramsel", tmp_paramsel);
        opt_tmp->addOpt("optimizer", tmp_optimizer);

        tmp_kernel->addOpt("type", "rbf");
        tmp_paramsel->addOpt("sigma", new OptNumber(sigma));


        PredKernelTrainTest<T> predKernelTask;
        GurlsOptionsList* predKernel = predKernelTask.execute(*Xva, empty, *opt_tmp);

        OptMatrix<gMat2D<T> > *k_opt = predKernel->getOptAs<OptMatrix<gMat2D<T> > >("K");
        k_opt->detachValue();
        predKernel_K = &(k_opt->getValue());

        delete predKernel;
        opt_tmp->removeOpt("optimizer", false);
        delete opt_tmp;
    }

    opt->addOpt("optimizer", tmp_optimizer);


    std::set<unsigned long> ireg;
//    for(unsigned long i=0; i<n_rank; ++i)
//        ireg.insert( static_cast<unsigned long>(gurls::round(std::pow(m, (i+1.0)/n_rank)))-1);

    unsigned long step = static_cast<unsigned long>(floor(m/n_rank));
    for(unsigned long i = m-1, count = 0; count < n_rank; i-=step, ++count)
        ireg.insert(i);


    GurlsOptionsList* perf_opt = new GurlsOptionsList("perf_opt");
    gMat2D<T>* alpha = new gMat2D<T>();
    T maxPerf = -std::numeric_limits<T>::max();
    unsigned long maxRank = 0;
    gMat2D<T> *perfs_mat = new gMat2D<T>(1, ireg.size());
    T *perfs = perfs_mat->getData();
    set(perfs, (T)-1, ireg.size());

    gMat2D<T> *times_mat = new gMat2D<T>(1, ireg.size());
    T *times = times_mat->getData();
    set(times, (T)0, ireg.size());


    gMat2D<T> *guesses_mat = new gMat2D<T>(1, ireg.size());
    T *guesses_mat_it = guesses_mat->getData();
    set(guesses_mat_it, (T)-1, ireg.size());


    T* G = new T[n*m];  //Cholesky factor
    set(G, (T)0.0, n*m);

    T* Q = new T[n*m];  //Q part of the QR decomposition
    set(Q, (T)0.0, n*m);

    T* RR = new T[m*m];  //inverse of R*R'
    set(RR, (T)0.0, m*m);

    T* yPvec = new T[n*t];
    copy(yPvec, y.getData(), n*t);


    boost::posix_time::ptime begin, end;
    boost::posix_time::time_duration diff;

    T prevPerf = 0;

    // ---- Cholesky decomposition ----
    T* diagG;
    T* newKcol;
    unsigned long* pVec;
    T normG;

    begin = boost::posix_time::microsec_clock::local_time();
    for(unsigned long i=0; i<m; ++i)
    {

        if(i==0)
        {
            //First iteration
            diagG = new T[n];
            set(diagG, (T)1.0, n);

            pVec = new unsigned long[n];
            for(unsigned long *it = pVec, *const end = pVec+n, i=0; it != end; ++it, ++i)
                *it = i;

            G[0] = 1.0;


        //     newKcol = exp(-1/sigma^2*square_distance(X(Pvec(2:n),:)',X(1,:)'));
            newKcol = computeNewKcol(X.getData(), n, d, sigma, pVec, 1, n);

        //    G(2:n,1) = newKcol;
            copy(G + 1, newKcol, n-1);

        //    diagG(2:n)=ones(n-1,1)-sum(G(2:n,1).^2,2);
            mult(newKcol, newKcol, newKcol, n-1);
            axpy(n-1, (T)-1.0, newKcol, 1, diagG+1, 1); // assumes diagG has been previously initialized to 1

        //    normG = norm(G(:,1));
            normG = sqrt(sumv(newKcol, n-1) + (G[0]*G[0]));
            delete [] newKcol;

        //    Q(:,1) = G(:,1) / normG;
            copy(Q, G, n);
            scal(n, (T)1.0/normG, Q, 1);

        //    RR(1,1) = 1/(normG^2);
            RR[0] = 1.0/ (normG*normG);
        }
        else // all other iterations
        {
            // find best new element
            unsigned long jast = std::max_element(diagG+i, diagG+n)-diagG;

            // updates permutation
    //        Pvec( [i jast] ) = Pvec( [jast i] );
            std::swap(*(pVec+i), *(pVec+jast));

            // swap y
            gurls::swap(t, yPvec+i, n, yPvec+jast, n);

            // updates all elements of G due to new permutation
    //        G([i jast],1:i)=G([ jast i],1:i);
            gurls::swap(i+1, G+i, n, G+jast, n);


            // updates all elements of Q due to new permutation
    //        Q([i jast],1:i-1) = Q([ jast i],1:i-1);
            gurls::swap(i, Q+i, n, Q+jast, n);


            // do the cholesky update
    //        G(i,i)=diagG(jast);
    //        G(i,i)=sqrt(G(i,i));
            const T G_ii = sqrt(diagG[jast]);
            G[(n*i)+i] = G_ii;


    //        newKcol = exp(-1/sigma^2*square_distance(X(Pvec((i+1):n),:)',X(Pvec(i),:)'));
            newKcol = computeNewKcol(X.getData(), n, d, sigma, pVec, i+1, n);

    //        G((i+1):n,i)=1/G(i,i)*( newKcol - G((i+1):n,1:(i-1))*(G(i,1:(i-1)))');
            const int rows = (n-(i+1));
            const unsigned long cols = i;

            T* G_ip1_n = new T[rows*(cols+1)]; // extracting +1 cols for future use

            for(T *G_it = G+(i+1), *Gi_it = G_ip1_n, *const Gi_end = Gi_it+rows*cols; Gi_it != Gi_end; G_it+=n, Gi_it+=rows)
                copy(Gi_it, G_it, rows);

            T* G_i= new T[cols];
            copy(G_i, G+i, cols, 1, n);

            const T beta = 1.0/G_ii;

            gemv(CblasNoTrans, rows, cols, -beta, G_ip1_n, rows, G_i, 1, beta, newKcol, 1);
            copy(G+(i+1)+(n*i), newKcol, rows);

            delete [] G_i;

            // updates diagonal elements
    //        diagG((i+1):n)=ones(n-i,1)-sum(G((i+1):n,1:i).^2,2  );
            copy(G_ip1_n+(rows*cols), newKcol, rows);
            mult(G_ip1_n, G_ip1_n, G_ip1_n, rows*(cols+1));

            delete [] newKcol;

            T* sums = new T[rows];
            sum_col(G_ip1_n, sums, rows, cols+1);
            set(diagG+(i+1), (T)1.0, rows);
            axpy(rows, (T)-1.0, sums, 1, diagG+(i+1), 1);

            delete [] sums;
            delete [] G_ip1_n;


            // performs QR decomposition
    //        Gcol = G(:,i);
            T* Gcol = new T[n];
            copy(Gcol, G+(n*i), n);

    //        Rcol = Q(:,1:(i-1))' * Gcol;
            T* Rcol = new T[cols];
            T* Q_sub = new T[n*cols];
            copy(Q_sub, Q, n*cols);
            gemv(CblasTrans, n, cols, (T)1.0, Q_sub, n, Gcol, 1, (T)0.0, Rcol, 1);


    //        Q(:,i) = Gcol - Q(:,1:(i-1)) * Rcol;
            T* Q_i = Q+(n*i);
            copy(Q_i, Gcol, n);
            gemv(CblasNoTrans, n, cols, (T)-1.0, Q_sub, n, Rcol, 1, (T)1.0, Q_i, 1);

            delete [] Gcol;
            delete [] Q_sub;

    //        Rii = norm(Q(:,i));
    //        Q(:,i) = Q(:,i) / Rii;
            normG = nrm2(n, Q_i, 1);
            scal(n, (T)1.0/normG, Q_i, 1);


            // updates
    //        RR(1:(i-1),i) = -(RR(1:(i-1),1:(i-1))*Rcol)./Rii;
            T* RR_sub = new T[cols*cols];
            for(T *R_it = RR, *Rs_it = RR_sub, *const R_end = R_it+(i*m); R_it != R_end; R_it+=m, Rs_it+=cols)
                copy(Rs_it, R_it, cols);

            T* RR_i = RR+(m*i);
            T* RR_Rcol = new T[cols];

            gemv(CblasNoTrans, cols, cols, (T)1.0, RR_sub, cols, Rcol, 1, (T)0.0, RR_Rcol, 1);
            copy(RR_i, RR_Rcol, cols);
            scal(cols, (T)-1.0/normG, RR_i, 1);

            delete[] RR_sub;

    //        RR(i,1:(i-1)) = RR(1:(i-1),i)';
            copy(RR+i, RR_i, cols, m, 1);

    //        RR(i,i) = (Rcol'*RR(1:(i-1),1:(i-1))*Rcol + 1)./(Rii^2);
            normG *= normG;
            RR[i+(i*m)] = (dot(cols, Rcol, 1, RR_Rcol, 1)+1) / normG;

            delete [] Rcol;
            delete [] RR_Rcol;
        }


        if(ireg.find(i) != ireg.end())
        {
            *guesses_mat_it++ = i+1;

//            vout.alpha = Q(:,1:i)*RR(1:i,1:i)*(Q(:,1:i)'*y(Pvec,:));

            const unsigned long ii = i+1;


            T* QtYp = new T[ii*t];
            dot(Q, yPvec, QtYp, n, ii, n, t, ii, t, CblasTrans, CblasNoTrans, CblasColMajor);

            T* RR_sub = new T[ii*ii];
            for(T *R_it = RR, *Rs_it = RR_sub, *const R_end = R_it+(ii*m); R_it != R_end; R_it+=m, Rs_it+=ii)
                copy(Rs_it, R_it, ii);

            T* RRQtYp = new T[ii*t];
            dot(RR_sub, QtYp, RRQtYp, ii, ii, ii, t, ii, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

            delete [] RR_sub;
            delete [] QtYp;

            T* alpha_i = new T[n*t];
            dot(Q, RRQtYp, alpha_i, n, ii, ii, t, n, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

            delete [] RRQtYp;



        //    permute indices of alpha
        //    vout.alpha(Pvec, :) = vout.alpha;
            gMat2D<T>* alphaMat = new gMat2D<T>(n,t);

            T *const alphaMat_it = alphaMat->getData();
            for (unsigned long j=0; j<n; ++j)
                copy(alphaMat_it+pVec[j], alpha_i+j, t, n, n);

            delete [] alpha_i;

            if(computePred)
            {
                // pred_primal
                gMat2D<T> * pred = new gMat2D<T>(Xva->rows(), t);
                dot(predKernel_K->getData(), alphaMat_it, pred->getData(), Xva->rows(), n, n, t, Xva->rows(), t, CblasNoTrans, CblasNoTrans, CblasColMajor);

                // perf_macroavg

                perf_opt->addOpt("pred", new OptMatrix<gMat2D<T> >(*pred));
                GurlsOptionsList* perf = perfTask->execute(empty, *yva, *perf_opt);

                gMat2D<T>& acc = perf->getOptValue<OptMatrix<gMat2D<T> > >("acc");
                T perf_i = sumv(acc.getData(), acc.getSize())/acc.getSize();
                *perfs = perf_i;
//                ++perfs;

                perf_opt->removeOpt("pred");
                delete perf;

                if(perf_i > maxPerf)
                {
                    maxPerf = perf_i;
                    maxRank = i;

                    delete alpha;
                    alpha = alphaMat;
                }
                else
                    delete alphaMat;


                if(le(*perfs, prevPerf) && i > 10)
                    break;
                else
                {
                    prevPerf = *perfs;
                    ++perfs;
                }
            }
            else
            {
                if(i == m-1)
                {
                    delete alpha;
                    alpha = alphaMat;
                    maxPerf = 1;
                    maxRank = i;

                    set(perfs, (T)0.0, ireg.size());
                    set(times, (T)0.0, ireg.size());
                }
            }

            end = boost::posix_time::microsec_clock::local_time();
            diff = end-begin;

            *times = diff.total_milliseconds();
            ++times;

        }

    }

    delete perf_opt;
    delete perfTask;

    delete [] G;
    delete [] diagG;


    delete [] yPvec;
    delete [] RR;
    delete [] Q;

    delete [] pVec;

    delete predKernel_K;

    GurlsOptionsList* paramsel = opt->getOptAs<GurlsOptionsList>("paramsel");
    paramsel->removeOpt("alpha");
    paramsel->removeOpt("acc");
    paramsel->removeOpt("maxRank");
    paramsel->removeOpt("maxPerf");
    paramsel->removeOpt("times");
    paramsel->removeOpt("guesses");

    paramsel->addOpt("alpha", new OptMatrix<gMat2D<T> >(*alpha));
    paramsel->addOpt("acc", new OptMatrix<gMat2D<T> >(*perfs_mat));
    paramsel->addOpt("maxRank", new OptNumber(maxRank));
    paramsel->addOpt("maxPerf", new OptNumber(maxPerf));
    paramsel->addOpt("times", new OptMatrix<gMat2D<T> >(*times_mat));
    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

}

template <typename T>
void ICholWrapper<T>::update(const gVec<T> &X, const gVec<T> &y)
{
}

template <typename T>
gMat2D<T>* ICholWrapper<T>::eval(const gMat2D<T> &X)
{
    GurlsOptionsList *opt = this->opt;

    const gMat2D<T> &alpha_mat = opt->getOptValue<OptMatrix<gMat2D<T> > >("paramsel.alpha");
    const T *const alpha = alpha_mat.getData();

    const unsigned long n = X.rows();
    const unsigned long nt = alpha_mat.rows();
    const unsigned long t = alpha_mat.cols();

    PredKernelTrainTest<T> predKernelTask;

    gMat2D<T> empty;


    GurlsOptionsList* opt_tmp = new GurlsOptionsList("tmp");
    GurlsOptionsList* tmp_kernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* tmp_paramsel = new GurlsOptionsList("paramsel");

    opt_tmp->addOpt("kernel", tmp_kernel);
    opt_tmp->addOpt("paramsel", tmp_paramsel);
    opt_tmp->addOpt("optimizer", opt->getOpt("optimizer"));

    tmp_kernel->addOpt("type", "rbf");
    tmp_paramsel->addOpt("sigma", new OptNumber(opt->getOptAsNumber("paramsel.sigma")));


    gMat2D<T> *y = new gMat2D<T>(n, t);

    GurlsOptionsList *predKernel = predKernelTask.execute(X, empty, *opt_tmp);
    const gMat2D<T> &predKernel_K = predKernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

    dot(predKernel_K.getData(), alpha, y->getData(), n, nt, nt, t, n, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

    opt_tmp->removeOpt("optimizer", false);
    delete opt_tmp;
    delete predKernel;

    return y;
}

template <typename T>
gMat2D<T>* ICholWrapper<T>::eval_ls(const gMat2D<T> &X)
{
    GurlsOptionsList *opt = this->opt;

    const gMat2D<T> &alpha_mat = opt->getOptValue<OptMatrix<gMat2D<T> > >("paramsel.alpha");
    const T *const alpha = alpha_mat.getData();

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long nt = alpha_mat.rows();
    const unsigned long t = alpha_mat.cols();

    gMat2D<T> Xi_mat(1, d);
    T* yi = new T[t];
    T* const Xi = Xi_mat.getData();

    gMat2D<T> *y_mat = new gMat2D<T>(n, t);
    T* y_it = y_mat->getData();

    GurlsOptionsList* opt_tmp = new GurlsOptionsList("tmp");
    GurlsOptionsList* tmp_kernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* tmp_paramsel = new GurlsOptionsList("paramsel");

    opt_tmp->addOpt("kernel", tmp_kernel);
    opt_tmp->addOpt("paramsel", tmp_paramsel);
    opt_tmp->addOpt("optimizer", opt->getOpt("optimizer"));

    tmp_kernel->addOpt("type", "rbf");
    tmp_paramsel->addOpt("sigma", new OptNumber(opt->getOptAsNumber("paramsel.sigma")));


    PredKernelTrainTest<T> predKernelTask;
    const gMat2D<T> empty;

    for(int i=0; i<n; ++i, ++y_it)
    {
        getRow(X.getData(), n, d, i, Xi);

        GurlsOptionsList *predKernel = predKernelTask.execute(Xi_mat, empty, *opt_tmp);
        const gMat2D<T> &predKernel_K = predKernel->getOptValue<OptMatrix<gMat2D<T> > >("K");

        dot(predKernel_K.getData(), alpha, yi, 1, nt, nt, t, 1, t, CblasNoTrans, CblasNoTrans, CblasColMajor);
        copy(y_it, yi, t, n, 1);

        delete predKernel;
    }

    opt_tmp->removeOpt("optimizer", false);
    delete opt_tmp;
    delete [] yi;

    return y_mat;
}

template <typename T>
void ICholWrapper<T>::setRankMax(unsigned long rank)
{
    GurlsOptionsList* paramsel = this->opt->template getOptAs<GurlsOptionsList>("paramsel");
    paramsel->removeOpt("rank_max");
    paramsel->addOpt("rank_max", new OptNumber(rank));
}

template <typename T>
void ICholWrapper<T>::setNRank(unsigned long n_rank)
{
    GurlsOptionsList* paramsel = this->opt->template getOptAs<GurlsOptionsList>("paramsel");
    paramsel->removeOpt("n_rank");
    paramsel->addOpt("n_rank", new OptNumber(n_rank));
}

template <typename T>
void ICholWrapper<T>::setSigma(double sigma)
{
    GurlsOptionsList* paramsel = this->opt->template getOptAs<GurlsOptionsList>("paramsel");
    paramsel->removeOpt("sigma");
    paramsel->addOpt("sigma", new OptNumber(sigma));
}

template <typename T>
void ICholWrapper<T>::setXva(const gMat2D<T> &Xva)
{
    GurlsOptionsList* split = this->opt->template getOptAs<GurlsOptionsList>("split");
    gMat2D<T> *Xva_mat = new gMat2D<T>(Xva);
    split->removeOpt("Xva");
    split->addOpt("Xva", new OptMatrix<gMat2D<T> >(*Xva_mat));
}

template <typename T>
void ICholWrapper<T>::setyva(const gMat2D<T> &yva)
{
    GurlsOptionsList* split = this->opt->template getOptAs<GurlsOptionsList>("split");
    gMat2D<T> *yva_mat = new gMat2D<T>(yva);
    split->removeOpt("yva");
    split->addOpt("yva", new OptMatrix<gMat2D<T> >(*yva_mat));
}

template <typename T>
T* ICholWrapper<T>::computeNewKcol(const T* Xtr, const unsigned long xr, const unsigned long xc,
                                             const double sigma,
                                             const unsigned long* pVec,
                                             const unsigned long start, const unsigned long n)
{
    // newKcol = exp(-1/sigma^2*square_distance(Xtr(Pvec(start:n),:)',Xtr(Pvec(start-1),:)'));

    const unsigned long len = n-start;

    T* newKcol = new T[len];

    T* XtrPvec = new T[len*xc];
    subMatrixFromRows(Xtr, xr, xc, pVec+start, len, XtrPvec);

    T* XtrPvecT = new T[xc*len];
    transpose(XtrPvec, len, xc, XtrPvecT);
    delete[] XtrPvec;

    T* XtrRow = new T[xc];
    getRow(Xtr, xr, xc, pVec[start-1], XtrRow);

    distance(XtrPvecT, XtrRow, xc, len, 1, newKcol);
    delete [] XtrRow;
    delete [] XtrPvecT;

    scal(len, (T)(-1.0/(sigma*sigma)), newKcol, 1);
    gurls::exp(newKcol, len);

    return newKcol;
}

}
