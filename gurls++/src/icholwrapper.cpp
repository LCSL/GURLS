#include "icholwrapper.h"

#include "optmatrix.h"
#include "predkerneltraintest.h"
#include "primal.h"
#include "utils.h"
#include "splitho.h"

using namespace gurls;

ICholWrapper::ICholWrapper(const std::string& name):GurlsWrapper(name)
{
    opt = new GurlsOptionsList(name, true);

    GurlsOptionsList *paramsel = new GurlsOptionsList("paramsel");
    opt->addOpt("paramsel", paramsel);
}

void ICholWrapper::train(const gMat2D<T> &X, const gMat2D<T> &y)
{
    const unsigned long m = static_cast<unsigned long>(opt->getOptAsNumber("paramsel.rank_opt"));
    const double sigma = opt->getOptAsNumber("paramsel.sigma");

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = y.cols();


    T* G = new T[n*m];  //Cholesky factor
    set(G, (T)0.0, n*m);

    T* Q = new T[n*m];  //Q part of the QR decomposition
    set(Q, (T)0.0, n*m);

    T* RR = new T[m*m];  //inverse of R*R'
    set(RR, (T)0.0, m*m);


    // ---- Cholesky decomposition ----

    //First iteration
    T* diagG = new T[n];
    set(diagG, (T)1.0, n);

    unsigned long* pVec = new unsigned long[n];
    for(unsigned long *it = pVec, *const end = pVec+n, i=0; it != end; ++it, ++i)
        *it = i;

    G[0] = 1.0;


//     newKcol = exp(-1/sigma^2*square_distance(X(Pvec(2:n),:)',X(1,:)'));
    T* newKcol = computeNewKcol(X.getData(), n, d, sigma, pVec, 1, n);

//    G(2:n,1) = newKcol;
    copy(G + 1, newKcol, n-1);

//    diagG(2:n)=ones(n-1,1)-sum(G(2:n,1).^2,2);
    mult(newKcol, newKcol, newKcol, n-1);
    axpy(n-1, (T)-1.0, newKcol, 1, diagG+1, 1); // assumes diagG has been previously initialized to 1

//    normG = norm(G(:,1));
    T normG = sqrt(sumv(newKcol, n-1) + (G[0]*G[0]));
    delete [] newKcol;

//    Q(:,1) = G(:,1) / normG;
    copy(Q, G, n);
    scal(n, (T)1.0/normG, Q, 1);

//    RR(1,1) = 1/(normG^2);
    RR[0] = 1.0/ (normG*normG);

    // all other iterations
    for(unsigned long i=1; i<m; ++i)
    {
        // find best new element
        unsigned long jast = std::max_element(diagG+i, diagG+n)-diagG;

        // updates permutation
//        Pvec( [i jast] ) = Pvec( [jast i] );
        std::swap(*(pVec+i), *(pVec+jast));


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


    delete [] G;
    delete [] diagG;


//    vout.alpha = Q*RR*(Q'*y(Pvec,:));

    T* yPvec = new T[n*t];
    subMatrixFromRows(y.getData(), n, t, pVec, n, yPvec);

    T* QtYp = new T[m*t];
    dot(Q, yPvec, QtYp, n, m, n, t, m, t, CblasTrans, CblasNoTrans, CblasColMajor);

    delete [] yPvec;

    T* RRQtYp = new T[m*t];
    dot(RR, QtYp, RRQtYp, m, m, m, t, m, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

    delete [] RR;
    delete [] QtYp;

    T* alpha = new T[n*t];
    dot(Q, RRQtYp, alpha, n, m, m, t, n, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

    delete [] Q;

    gMat2D<T>* alphaMat = new gMat2D<T>(n,t);

//    permute indices of alpha
//    vout.alpha(Pvec, :) = vout.alpha;

    T*const alphaMat_it = alphaMat->getData();
    for (unsigned long i=0; i<n; ++i)
        copy(alphaMat_it+pVec[i], alpha+i, t, n, n);

    delete [] pVec;
    delete [] alpha;


    GurlsOptionsList* paramsel = opt->getOptAs<GurlsOptionsList>("paramsel");
    paramsel->addOpt("alpha", new OptMatrix<gMat2D<T> >(*alphaMat));
}

void ICholWrapper::update(const gVec<T> &X, const gVec<T> &y)
{
}

gMat2D<GurlsWrapper::T>* ICholWrapper::eval(const gMat2D<T> &X)
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

void ICholWrapper::setRank(unsigned long rank)
{
    GurlsOptionsList* paramsel = opt->getOptAs<GurlsOptionsList>("paramsel");
    paramsel->addOpt("rank_opt", new OptNumber(rank));
}

void ICholWrapper::setSigma(double sigma)
{
    GurlsOptionsList* paramsel = opt->getOptAs<GurlsOptionsList>("paramsel");
    paramsel->addOpt("sigma", new OptNumber(sigma));
}


GurlsWrapper::T* ICholWrapper::computeNewKcol(const T* Xtr, const unsigned long xr, const unsigned long xc,
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
