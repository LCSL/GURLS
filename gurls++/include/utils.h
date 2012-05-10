/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, IIT@MIT Lab
 * All rights reserved.
 *
 * author:  M. Santoro
 * email:   msantoro@mit.edu
 * website: http://cbcl.mit.edu/IIT@MIT/IIT@MIT.html
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name(s) of the copyright holders nor the names
 *       of its contributors or of the Massacusetts Institute of
 *       Technology or of the Italian Institute of Technology may be
 *       used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef _GURLS_UTILS_H_
#define _GURLS_UTILS_H_

#include <map>
#include <algorithm>

#include "gmath.h"
#include "gmat2d.h"

#include "optlist.h"
#include "primal.h"
#include "macroavg.h"

using namespace std;

#define FLT_EQUALS(flt1, flt2) \
        ( flt1-flt2 >= -FLT_EPSILON && flt1-flt2 <= FLT_EPSILON )

#define FLT_GT(a, b)\
    ((a - b) > ( std::min(fabs(a), fabs(b))* FLT_EPSILON))

#define GT(a, b)\
    ((a - b) > ( std::min(fabs(a), fabs(b))* std::numeric_limits<T>::epsilon()))

namespace gurls {

/*
 * compute the average precision through precision and recall
 */
template <typename T>
double precrec_driver(T* out, T* gt, unsigned long N)
{

    multimap<T, T> data;
    for (unsigned long i = 0; i < N; i++){
        data.insert(pair<T,T>(-out[i], gt[i]));
    }

    T tp[N], fp[N], prec[N], rec[N];
    unsigned long tpcumsum = 0;
    unsigned long fpcumsum = 0;

    typename multimap<T, T>::iterator it = data.begin();
    typename multimap<T, T>::iterator end = data.end();


    unsigned long idx = 0;
    while (it != end) {
        tp[idx] = it->second > 0 ? tpcumsum+=1 : tpcumsum;
        fp[idx] = it->second < 0 ? fpcumsum+=1 : fpcumsum;
        it++;
        idx++;
    }

    for (idx = 0; idx < N; idx++){
        rec[idx] = tp[idx]/tpcumsum;
        prec[idx] = tp[idx]/(tp[idx]+fp[idx]);
    }

    // compute average precision
    T ap = 0;
    T p = 0.0;
    T t = 0.0;

    const int stepsNumber = 11;

    for(int steps = 0; steps<stepsNumber; ++steps)
    {
        p = 0.0;
        for (idx = 0; idx < N; idx++){
            if (FLT_GT(rec[idx], t) || FLT_EQUALS(rec[idx], t)){
                p = std::max(p, prec[idx]);
            }
        }
        ap+=p;
        t+=.1;
    }

    return ap/=stepsNumber;
}


template <typename T>
void rls_pegasos_driver(const T* X, const T* bY, GurlsOptionsList& opt,
                        const int X_rows, const int X_cols,
                        const int bY_rows, const int bY_cols)
{
    //	lambda = opt.singlelambda(opt.paramsel.lambdas);
    std::vector<double> ll = OptNumberList::dynacast(opt.getOpt("lambdas"))->getValue();
    T lambda = static_cast<T>((OptFunction::dynacast(opt.getOpt("singlelambda")))->getValue(ll.data(), ll.size()));

//            %% Inputs

//            [n,d] = size(X);
    const int n = X_rows;
    const int d = X_cols;

//            [T] = size(bY,2);
    const int t = bY_cols;

//            %% Initialization
//            cfr = opt.cfr;

//            W = cfr.W; %dxT
    GurlsOption *W_opt = opt.getOpt("W");
    gMat2D<T> *W_mat = &(OptMatrix<gMat2D<T> >::dynacast(W_opt))->getValue();
    gMat2D<T> W(W_mat->cols(), W_mat->rows());
    W_mat->transpose(W);
//    T* W = new T[d*t];

//            W_sum = cfr.W_sum;
    GurlsOption *W_sum_opt = opt.getOpt("W_sum");
    gMat2D<T> *W_sum_mat = &(OptMatrix<gMat2D<T> >::dynacast(W_sum_opt))->getValue();
    gMat2D<T> W_sum(W_sum_mat->cols(), W_sum_mat->rows());
    W_sum_mat->transpose(W_sum);

//            count = cfr.count; %contatore per numero iterazioni totali, prende quello in input e lo aggiorna
    int count = static_cast<int>(opt.getOptAsNumber("count"));

//            t0 = cfr.t0;
    T t0 = static_cast<T>(opt.getOptAsNumber("t0"));


    unsigned int * seq = new unsigned int[n];
    T* xt = new T[d];
    T* y_hat = new T[t];
    T* r = new T[t];
    const int W_size = d*t;
    T* xtr = new T[W_size];

//            %% Initialization
//            iter = 0;
    int iter;

    const T thr = sqrt(t/lambda);

//            seq = randperm(n); %riordina 1:n in modo casuale
    randperm(n, seq);

    for(iter = 0; iter<n; ++iter)
    {
//            while iter < n,
//                iter = iter + 1;
//                idx = seq(iter);
        const unsigned int idx = seq[iter]-1;

//                %% Stepsize

//                %% Update Equations
//                xt = X(idx,:); %1xd
        getRow(X, n, d, idx, xt);

//                y_hat = (xt*W); %1xT
        dot(xt, W.getData(), y_hat, 1, d, d, t, 1, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

//                r = bY(idx,:) - y_hat; %1xT
        getRow(bY, bY_rows, t, idx, r);
        axpy(t, (T)-1.0, y_hat, 1, r, 1);


//                eta = 1.0/(lambda*(count + t0));
        const T eta = ((T)1.0)/(lambda*(count + t0));

//                W = (1 - lambda*eta)*W + eta*xt'*r; %dxT
        const T coeff = 1.0 - (lambda*eta);
//        dot(xt, r, xtr, d, 1, 1, t, d, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

//        set(xtr, (T)0.0, W_size);
//        axpy(W_size, coeff, W.getData(), 1, xtr, 1);
//        copy(W.getData(), xtr, W_size);
        scal(W_size, coeff, W.getData(), 1);

        dot(xt, r, xtr, d, 1, 1, t, d, t, CblasNoTrans, CblasNoTrans, CblasColMajor);
        axpy(W_size, eta, xtr, 1, W.getData(), 1);


//                %% Projection onto the ball with radius sqrt(T/lambda)
//                nW = norm(W,'fro'); %Frobenius norm

        float tmp = dot(W_size, W.getData(), 1, W.getData(), 1);
        T  nW = sqrt(tmp);


//                if nW > sqrt(T/lambda)
        if( GT(nW, thr) )
        {
//                    W = (W/nW)*sqrt(T/lambda);
            set(xtr, (T)0.0, W_size);
            axpy(W_size, thr/nW, W.getData(), 1, xtr, 1);
            copy(W.getData(), xtr, W_size);
        }

//                %% Averaging

//                W_sum = W_sum + W;
        axpy(W_size, (T)1.0, W.getData(), 1, W_sum.getData(), 1);

//                count = count + 1;
        ++count;

//                %% Testing
//                if(mod(count,n) == 1)
//                    fprintf('\n\tObjective : %f',obj_primal(W, X, bY, lambda));
//                    cfr.acc_last(end+1) = test_classifier (W,opt);
//                    fprintf('\n\tLast Acc: %f', cfr.acc_last(end));
//                    cfr.acc_avg(end+1) = test_classifier (W_sum/count,opt);
//                    fprintf('\n\tAvg Acc: %f\n', cfr.acc_avg(end));
//                end
    }

    delete[] seq;
    delete[] xt;
    delete[] y_hat;
    delete[] r;
    delete[] xtr;

//            cfr.W = W;
    W.transpose(*W_mat);

//            cfr.W_last = W;
    if(opt.hasOpt("W_last"))
    {
        GurlsOption *W_last_opt = opt.getOpt("W_last");
        gMat2D<T> *W_last = &(OptMatrix<gMat2D<T> >::dynacast(W_last_opt))->getValue();
        copy(W_last->getData(), W_mat->getData(), W_size);
    }
    else
    {
        gMat2D<T> *W_last = new gMat2D<T>(*W_mat);
        opt.addOpt("W_last", new OptMatrix<gMat2D<T> >(*W_last));
    }

//            cfr.W_sum = W_sum;
    W_sum.transpose(*W_sum_mat);

//            cfr.count = count;
    opt.removeOpt("count");
    opt.addOpt("count", new OptNumber(count));

//            cfr.iter = iter;
    opt.removeOpt("iter");
    opt.addOpt("iter", new OptNumber(iter));

    //	cfr.C = [];
    opt.removeOpt("C");
    gMat2D<T>* emptyC = new gMat2D<T>();
    opt.addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

    //	cfr.X = [];
    opt.removeOpt("X");
    gMat2D<T>* emptyX = new gMat2D<T>();
    opt.addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));

}

// function [acc] = test_classifier (W, opt)
template <typename T>
T test_classifier(const T* W, GurlsOptionsList& opt, const int rows, const int cols)
{
    //opt.rls.W = W;
    //copy W
    //copy<>()
//     GurlsOptionsList* rls = static_cast<GurlsOptionsList*>(opt.getOpt("rls"));
//     if(rls->hasOpt("W"))
//       rls->removeOpt("W");
//     rls->addOpt("W", new OptMatrix< gMat2D<T> > (*W));
    //opt.addOpt("rls",rls);
    gMat2D<T> *W_mat = NULL;
    if(!opt.hasOpt("W"))
    {
      W_mat = new gMat2D<T>(rows,cols);
      opt.addOpt("W", new OptMatrix<gMat2D<T> >(*W_mat));
    }

    if(W_mat == NULL)
    {
     GurlsOption *W_opt = opt.getOpt("W");
     W_mat = &(OptMatrix<gMat2D<T> >::dynacast(W_opt))->getValue();
    }

    gMat2D<T> W_t(W,W_mat->cols(),W_mat->rows(),false);
    W_t.transpose(*W_mat);

    GurlsOption *x = opt.getOpt("Xte");
    gMat2D<T>* X = &(OptMatrix< gMat2D<T> >::dynacast(x))->getValue();
    GurlsOption *y = opt.getOpt("yte");
    gMat2D<T>* Y = &(OptMatrix< gMat2D<T> >::dynacast(y))->getValue();

    //opt.pred = pred_primal(opt.Xte, opt.yte, opt);
    PredPrimal<T> pp;
    pp.execute(*X, *Y, opt);
    //opt.perf   = perf_macroavg(opt.Xte, opt.yte, opt);
    MacroAvg<T> ma;
    ma.execute(*X, *Y, opt);

    //acc = mean([opt.perf.acc]);
    GurlsOptionsList* perf = static_cast<GurlsOptionsList*>(opt.getOpt("perf"));
    GurlsOption *perf_opt = perf->getOpt("acc");
    gMat2D<T> *acc_mat = &(OptMatrix<gMat2D<T> >::dynacast(perf_opt))->getValue();
    const long size = acc_mat->getSize();
    T res;
    mean<T>(acc_mat.getData(),&res,1,size,1);
    return res;
}

}


#endif // _GURLS_UTILS_H_
