 /*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011, IIT@MIT Lab
  * All rights reserved.
  *
  * authors:  M. Santoro
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


#ifndef _GURLS_LOOCVPRIMAL_H_
#define _GURLS_LOOCVPRIMAL_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <set>

#include "options.h"
#include "optlist.h"
#include "gmat2d.h"
#include "gvec.h"
#include "gmath.h"

#include "paramsel.h"

#include "precisionrecall.h"
#include "macroavg.h"

namespace gurls {

template <typename T>
class LoocvPrimal: public ParamSelection<T>{

public:
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
    void executeLB_old(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
    void executeLB(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void LoocvPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt){

    executeLB(X, Y, opt);
    return;

    try {

        //[n,T]  = size(y);
        T n = Y.rows();
        T t = Y.cols();

        // This will be used in the following as temporary matrix
        gMat2D<T> *tmp, *tmp1;
        gVec<T> *tmpvec  = new gVec<T>(X.cols());
        *tmpvec = 0;

        // K = X'*X;
        gMat2D<T> Xt(X.cols(), X.rows());
        Xt = 0;
        X.transpose(Xt);
        gMat2D<T> K(Xt.rows(), X.cols());
        dot(Xt, X, K);

        // tot = opt.nlambda;
        T tot = std::ceil( static_cast<OptNumber*>(opt.getOpt("nlambda"))->getValue() );

        //	[Q,L] = eig(K);
        //	L = diag(L);
        gMat2D<T> Q(K.rows(), K.cols());
        gVec<T> L(K.cols());


//		std::cout << "K = "<< std::endl << K << std::endl;
        eig(K, Q, L);

//		std::cout << "K = "<< std::endl << K << std::endl;
//		std::cout << "Q = " << std::endl << Q << std::endl;

//		for (int i = 0; i < K.cols(); i++) {
//			std::cout << " " << L[i];
//		}
//		std::cout << std::endl;
        //gMat2D<T> Lmat = gMat2D<T>::diag(L);

        // WARNING ============================== -> SEE LAMBDAGUESSES
        //	filtered = L(L > 200*eps^0.5);
        // WARNING: using an approximate version of the eps Matlab variable
        //gVec<T> filtered = (L.asMatrix()).where(L.asMatrix() > 2.9802e-6f);
        gVec<T> filtered = L;
        // COMING SOON: a specialized version of the logical operators and of
        // the `where' method will be available for the gVec class very soon.

        //	lmin = min(filtered)/n;
        //	lmax = max(filtered)/n;
        T lmin = filtered.min()/ n;
        T lmax = filtered.max()/ n;

        //	q = (lmax/lmin)^(1/tot);
        T q = pow((lmax/lmin), ( static_cast<T>(1.0) / tot));

        //	guesses = zeros(1,tot);
        //	LOOSQE = zeros(tot,T);
        gVec<T> guesses(tot);
        guesses = 0;
        gMat2D<T> *LOOSQE = new gMat2D<T>(tot, t);
        *LOOSQE = 0;

        //	LEFT = X*Q;
        //	RIGHT = Q'*X'*y;
        gMat2D<T> LEFT(X.rows(), Q.cols());
        dot(X,Q,LEFT);
        gMat2D<T> Qt(Q.cols(), Q.rows());
        Q.transpose(Qt);
        tmp = new gMat2D<T> (Xt.rows(), Y.cols());
        dot(Xt, Y, *tmp);
        gMat2D<T> RIGHT(Qt.rows(), tmp->cols());
        dot(Qt, *tmp, RIGHT);
        delete tmp;

        //	right = Q'*X';
        gMat2D<T> right(Qt.rows(), Xt.cols());
        dot(Qt, Xt, right);
        // gMat2D<T> right(X*Q);


        //gMat2D<T> num(n,t);
        gVec<T>	  den(n);
        gMat2D<T> Le(n, t);

        //	for i = 1:tot
        for(int s = 0; s < tot; s++)
        {
            //	guesses(i) = lmin*(q^i);
            guesses[s] = lmin*std::pow(q,s+1);

            //		LL = L + (n*guesses(i));
            //		LL = LL.^(-1)
            //		LL = diag(LL);
            tmpvec = new gVec<T>(L+(n*guesses[s]));

            //*tmpvec = L+(n*guesses[s]);
            gMat2D<T> LL = gMat2D<T>::diag( tmpvec->reciprocal() );

            // WARNING: salvare tempo e memoria evitando di inizializzare
            // ad ogni passo
            tmp = new gMat2D<T>(LL.rows(), RIGHT.cols());
            tmp1 = new gMat2D<T>(LEFT.rows(), tmp->cols());

            dot(LL, RIGHT, *tmp);
            dot(LEFT, *tmp, *tmp1);

            //		num = y - LEFT*LL*RIGHT;
            gMat2D<T> num = (Y - *tmp1);

            delete tmp, tmp1;

            // den = zeros(n,1);
            den = 0;

            //		for j = 1:n
            //			den(j) = 1-LEFT(j,:)*LL*right(:,j);
            //		end

            tmp = new gMat2D<T>(LL.rows(), 1);
//			std::cout << std::endl << "got here!" << std::endl;

//			std::cout << LEFT.rows() << " - " << LEFT.cols() << std::endl;
//			std::cout << LL.rows() << " - " << LL.cols() << std::endl;
//			std::cout << tmp->rows() << " - " << tmp->cols() << std::endl;

            for (int j = 0; j < n; j++){
                dot(LL, right(j).asMatrix(), *tmp);
                den[j] =  1 - dot(LEFT[j],(*tmp)(0));
            }

            //		Le = zeros(n,T);
            Le = 0;

            //		for t = 1:T
            //			Le(:,t) = num(:,t)./den;
            //		end
            for (unsigned long int k = 0; k < t; k++){
                tmpvec = new gVec<T>(num(k)*(den.reciprocal()));
                Le.setColumn( *tmpvec, k);
            }

            //		LOOSQE(i,:) = sum(Le.*Le);
            LOOSQE->setRow((Le*Le).sum(gurls::COLUMNWISE),s);

        }
        delete tmpvec;

        //	[dummy,bL] = min(LOOSQE);
        gVec<T>* bL = &LOOSQE->argmin(gurls::COLUMNWISE);
        //	vout.lambdas = guesses(bL);
        gVec<T>* lambdas = &guesses.copyLocations(*bL);

        OptNumberList* LAMBDA = new OptNumberList();
        for (int i = 0; i < lambdas->getSize(); i++) {
            LAMBDA->add(static_cast<double>(lambdas->at(i)));
        }
        opt.addOpt("lambdas", LAMBDA);
        //	vout.guesses = guesses;
        //opt.addOpt("guesses", new OptMatrix<T>(guesses.asMatrix()));
        opt.addOpt("guesses", new OptMatrix<gMat2D<T> >(guesses.asMatrix()));
        //	vout.LOOSQE = LOOSQE;
        //opt.addOpt("LOOSQE", new OptMatrix<T>(*LOOSQE));
        opt.addOpt("LOOSQE", new OptMatrix<gMat2D<T> >(*LOOSQE));
//		std::cout << "LOOSQE = " << std::endl << LOOSQE << std::endl;
//		std::cout << "guesses =" << std::endl << guesses << std::endl;
//		std::cout << "lambdas = " << std::endl << *lambdas << std::endl;

    }catch( gException& e){
        throw e;
    }

    return;
}

template <typename T>
void LoocvPrimal<T>::executeLB_old(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt){

    typename std::set<T*> garbage;

    try
    {

        //[n,T]  = size(y);
        const unsigned long n = Y_OMR.rows();
        const unsigned long t = Y_OMR.cols();

        gMat2D<float> X(X_OMR.cols(), X_OMR.rows());
        X_OMR.transpose(X);

        gMat2D<float> Y(Y_OMR.cols(), Y_OMR.rows());
        Y_OMR.transpose(Y);


        //	K = X'*X;
        const unsigned long xc = X_OMR.cols();
        const unsigned long xr = X_OMR.rows();

        if(xr != n)
            throw new gException("X and Y must have the same row count");


        T* K = new T[xc*xc];
        garbage.insert(K);
        dot(X.getData(), X.getData(), K, xr, xc, xr, xc, xc, xc, CblasTrans, CblasNoTrans, CblasColMajor);

        //		std::cout << "K = "<< std::endl << K << std::endl;


        // tot = opt.nlambda;
        int tot = std::ceil( static_cast<OptNumber*>(opt.getOpt("nlambda"))->getValue() );

        //	[Q,L] = eig(K);
        //	L = diag(L);

        T* Q, *L;

        eig(K, Q, L, xc, xc);
        garbage.insert(Q);
        garbage.insert(L);


//		std::cout << "K = "<< std::endl << K << std::endl;
//		std::cout << "Q = " << std::endl << Q << std::endl;

//		for (int i = 0; i < K.cols(); i++) {
//			std::cout << " " << L[i];
//		}
//		std::cout << std::endl;
        delete[] K;
        garbage.erase(K);
        //gMat2D<T> Lmat = gMat2D<T>::diag(L);

        // WARNING ============================== -> SEE LAMBDAGUESSES
        //	filtered = L(L > 200*eps^0.5);
        // WARNING: using an approximate version of the eps Matlab variable
        //gVec<T> filtered = (L.asMatrix()).where(L.asMatrix() > 2.9802e-6f);
//        gVec<T> filtered = L;
        T* filtered = L;
        // COMING SOON: a specialized version of the logical operators and of
        // the `where' method will be available for the gVec class very soon.

        //	lmin = min(filtered)/n;
        //	lmax = max(filtered)/n;
//        const T lmin = (*std::min_element(filtered, filtered+xc))/ n;
//        const T lmax = (*std::max_element(filtered, filtered+xc))/ n;

        //	q = (lmax/lmin)^(1/tot);
//        const T q = pow((lmax/lmin), ( static_cast<T>(1.0) / tot));

        //	guesses = zeros(1,tot);
        //	LOOSQE = zeros(tot,T);
//        T* guesses = new T[tot];
        T* guesses = lambdaguesses(filtered, xc, xc, xr, tot, (T)(opt.getOptAsNumber("smallnumber")));
        garbage.insert(guesses);
//        set(guesses, 0.f, tot);

        T* LOOSQE = new T[tot*t];
        garbage.insert(LOOSQE);
        set(LOOSQE, 0.f, tot*t);

        //	LEFT = X*Q;
        //	RIGHT = Q'*X'*y;
        T* LEFT = new T[xr*xc];
        garbage.insert(LEFT);
        dot(X.getData(), Q, LEFT, xr, xc, xc, xc, xr, xc, CblasNoTrans, CblasNoTrans, CblasColMajor);

        T* tmp = new T[xc*t];
        garbage.insert(tmp);
        dot(X.getData(), Y.getData(), tmp, xr, xc, n, t, xc, t, CblasTrans, CblasNoTrans, CblasColMajor);

        T* RIGHT = new T[xc*t];
        garbage.insert(RIGHT);
        dot(Q, tmp, RIGHT, xc, xc, xc, t, xc, t, CblasTrans, CblasNoTrans, CblasColMajor);

        delete[] tmp;
        garbage.erase(tmp);

        //	right = Q'*X';
        T* right = new T[xc*xr];
        garbage.insert(right);
        dot(Q, X.getData(), right, xc, xc, xr, xc, xc, xr, CblasTrans, CblasTrans, CblasColMajor);

        delete[] Q;
        garbage.erase(Q);

        T* den = new T[n];
        T* Le = new T[n*t];
        garbage.insert(den);
        garbage.insert(Le);

        T* tmpvec, *tmp1;

        //	for i = 1:tot
        for(int s = 0; s < tot; ++s)
        {
            //	guesses(i) = lmin*(q^i);
//            guesses[s] = lmin*std::pow(q,s+1);

            //		LL = L + (n*guesses(i));
            tmpvec = new T[xc];
            garbage.insert(tmpvec);
            set(tmpvec, n*guesses[s] , xc);
            axpy(xc, 1.0f, L, 1, tmpvec, 1);

            //		LL = LL.^(-1)
            setReciprocal(tmpvec, xc);
            //		LL = diag(LL);
            T* LL = diag(tmpvec, xc);
            garbage.insert(LL);

            delete[] tmpvec;
            garbage.erase(tmpvec);

            // WARNING: salvare tempo e memoria evitando di inizializzare
            // ad ogni passo

            //		num = y - LEFT*LL*RIGHT;
            tmp = new T[xc*t];
            garbage.insert(tmp);
            dot(LL, RIGHT, tmp, xc, xc, xc, t, xc, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

            tmp1 = new T[xr*t];
            garbage.insert(tmp1);
            dot(LEFT, tmp, tmp1, xr, xc, xc, t, xr, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

            delete[] tmp;
            garbage.erase(tmp);

            T* num = new T[xr*t];
            garbage.insert(num);
            // ?????
//            copy(num, Y.getData(), xr*t);
//            axpy(xr*t, -1.0f, tmp1, 1, num, 1);

            copy(num, tmp1, xr*t);
            axpy(xr*t, -1.0f, Y.getData(), 1, num, 1);

            delete[] tmp1;
            garbage.erase(tmp1);

            // den = zeros(n,1);
            set(den, 0.0f, n);

            //		for j = 1:n
            //			den(j) = 1-LEFT(j,:)*LL*right(:,j);
            //		end

            tmp = new T[xc];
            garbage.insert(tmp);
//			std::cout << std::endl << "got here!" << std::endl;

//			std::cout << LEFT.rows() << " - " << LEFT.cols() << std::endl;
//			std::cout << LL.rows() << " - " << LL.cols() << std::endl;
//			std::cout << tmp->rows() << " - " << tmp->cols() << std::endl;

            T* row = new T[xc];
            garbage.insert(row);

            for (int j = 0; j < n; ++j)
            {
                dot(LL, right +(xc*j), tmp, xc, xc, xc, 1, xc, 1, CblasNoTrans, CblasNoTrans, CblasColMajor);

                //extract j-th row from LEFT
                copy(row, LEFT + j, xc, 1, xr);

                den[j] =  ((T) 1.0) - dot (xc, row, 1, tmp, 1);
            }

            delete[] row;
            garbage.erase(row);
            delete[] tmp;
            garbage.erase(tmp);
            delete[] LL;
            garbage.erase(LL);

            //		Le = zeros(n,T);
            //set(Le, 0.0f, n*t);

            //		for t = 1:T
            //			Le(:,t) = num(:,t)./den;
            //		end
            setReciprocal(den, n);
            for (unsigned long int k = 0; k < t; ++k)
                mult(num + (xr*k), den, Le + (n*k), xr);

            //		LOOSQE(i,:) = sum(Le.*Le);
            mult(Le, Le, Le , n*t);
            row = new T[t];
            garbage.insert(row);
            sum(Le, row, n, t, t);

            copy(LOOSQE + s, row, t, tot, 1);

            delete[] row;
            garbage.erase(row);

            delete[] num;
            garbage.erase(num);
        }

        delete[] L;
        garbage.erase(L);
        delete[] LEFT;
        garbage.erase(LEFT);
        delete[] RIGHT;
        garbage.erase(RIGHT);
        delete[] right;
        garbage.erase(right);
        delete[] den;
        garbage.erase(den);
        delete[] Le;
        garbage.erase(Le);

        //	[dummy,bL] = min(LOOSQE);
        unsigned int* bL = new unsigned int[t];
        garbage.insert((T*)bL);
        argmin(LOOSQE, bL, tot, t, t);

        //	vout.lambdas = guesses(bL);
        T* lambdas = copyLocations(bL, guesses, t, tot);
        garbage.insert(lambdas);

        delete[] bL;
        garbage.erase((T*)bL);

        OptNumberList* LAMBDA = new OptNumberList();
        for (T* l_it = lambdas, *l_end = lambdas+t; l_it != l_end; ++l_it)
            LAMBDA->add(static_cast<double>(*l_it));

        delete[] lambdas;
        garbage.erase(lambdas);

        opt.addOpt("lambdas", LAMBDA);

        //	vout.guesses = guesses;
        //opt.addOpt("guesses", new OptMatrix<T>(guesses.asMatrix()));
        gVec<T> guessesVector(guesses, tot, false);
        opt.addOpt("guesses", new OptMatrix<gMat2D<T> >(guessesVector.asMatrix()));

        delete[] guesses;
        garbage.erase(guesses);

        //	vout.LOOSQE = LOOSQE;
        //opt.addOpt("LOOSQE", new OptMatrix<T>(*LOOSQE));
        tmp = transpose_rm(LOOSQE, t, tot);
        garbage.insert(tmp);

        delete[] LOOSQE;
        garbage.erase(LOOSQE);

        gMat2D<T>* LOOSQE_MAT = new gMat2D<T>(tmp, tot, t, true);

        delete[] tmp;
        garbage.erase(tmp);

        opt.addOpt("LOOSQE", new OptMatrix<gMat2D<T> >(*LOOSQE_MAT));
//		std::cout << "LOOSQE = " << std::endl << LOOSQE << std::endl;
//		std::cout << "guesses =" << std::endl << guesses << std::endl;
//		std::cout << "lambdas = " << std::endl << *lambdas << std::endl;

    }
    catch( gException& e)
    {

        for(typename std::set<T*>::iterator it = garbage.begin(); it != garbage.end(); ++it)
            delete[] (*it);

        throw e;
    }
}

template <typename T>
void LoocvPrimal<T>::executeLB(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt){

    typename std::set<T*> garbage;

    try
    {

        //[n,T]  = size(y);
        const unsigned long n = Y_OMR.rows();
        const unsigned long t = Y_OMR.cols();

        gMat2D<float> X(X_OMR.cols(), X_OMR.rows());
        X_OMR.transpose(X);

        gMat2D<float> Y(Y_OMR.cols(), Y_OMR.rows());
        Y_OMR.transpose(Y);


        //	K = X'*X;
        const unsigned long xc = X_OMR.cols();
        const unsigned long xr = X_OMR.rows();

        if(xr != n)
            throw new gException("X and Y must have the same row count");


        T* K = new T[xc*xc];
        garbage.insert(K);
        dot(X.getData(), X.getData(), K, xr, xc, xr, xc, xc, xc, CblasTrans, CblasNoTrans, CblasColMajor);

        //		std::cout << "K = "<< std::endl << K << std::endl;


        // tot = opt.nlambda;
        int tot = std::ceil( static_cast<OptNumber*>(opt.getOpt("nlambda"))->getValue() );

        //	[Q,L] = eig(K);
        //	L = diag(L);

//        T* Q, *L;

//        eig(K, Q, L, xc, xc);
//        garbage.insert(Q);
//        garbage.insert(L);


        T* Q = K;
        T* L = new T[xc];
        garbage.insert(L);

        eig_sm(Q, L, xc, xc);


//        delete[] K;
//        garbage.erase(K);

        // WARNING ============================== -> SEE LAMBDAGUESSES
        //	filtered = L(L > 200*eps^0.5);
        // WARNING: using an approximate version of the eps Matlab variable
        //gVec<T> filtered = (L.asMatrix()).where(L.asMatrix() > 2.9802e-6f);
//        gVec<T> filtered = L;
        T* filtered = L;
        // COMING SOON: a specialized version of the logical operators and of
        // the `where' method will be available for the gVec class very soon.

        //	lmin = min(filtered)/n;
        //	lmax = max(filtered)/n;
//        const T lmin = (*std::min_element(filtered, filtered+xc))/ n;
//        const T lmax = (*std::max_element(filtered, filtered+xc))/ n;

        //	q = (lmax/lmin)^(1/tot);
//        const T q = pow((lmax/lmin), ( static_cast<T>(1.0) / tot));

        //	guesses = zeros(1,tot);
        //	LOOSQE = zeros(tot,T);
//        T* guesses = new T[tot];
        T* guesses = lambdaguesses(filtered, xc, xc, xr, tot, (T)(opt.getOptAsNumber("smallnumber")));
        garbage.insert(guesses);
//        set(guesses, 0.f, tot);

        T* LOOSQE = new T[tot*t];
        garbage.insert(LOOSQE);
        set(LOOSQE, 0.f, tot*t);

        //	LEFT = X*Q;
        //	RIGHT = Q'*X'*y;
        T* LEFT = new T[xr*xc];
        garbage.insert(LEFT);
        dot(X.getData(), Q, LEFT, xr, xc, xc, xc, xr, xc, CblasNoTrans, CblasNoTrans, CblasColMajor);

        T* tmp = new T[xc*t];
        garbage.insert(tmp);
        dot(X.getData(), Y.getData(), tmp, xr, xc, n, t, xc, t, CblasTrans, CblasNoTrans, CblasColMajor);

        T* RIGHT = new T[xc*t];
        garbage.insert(RIGHT);
        dot(Q, tmp, RIGHT, xc, xc, xc, t, xc, t, CblasTrans, CblasNoTrans, CblasColMajor);

        delete[] tmp;
        garbage.erase(tmp);

        //	right = Q'*X';
        T* right = new T[xc*xr];
        garbage.insert(right);
        dot(Q, X.getData(), right, xc, xc, xr, xc, xc, xr, CblasTrans, CblasTrans, CblasColMajor);

        delete[] Q;
        garbage.erase(Q);

        T* den = new T[n];
//        T* Le = new T[n*t];
        garbage.insert(den);
//        garbage.insert(Le);

        T* tmpvec, *tmp1;

        GurlsOption* pred_old = NULL;

        if(opt.hasOpt("pred"))
        {
            pred_old = opt.getOpt("pred");
            opt.removeOpt("pred", false);
        }

        gMat2D<T>* pred = new gMat2D<T>(n, t);
        OptMatrix<gMat2D<T> >* pred_opt = new OptMatrix<gMat2D<T> >(*pred);
        opt.addOpt("pred", pred_opt);

        GurlsOptionsList* perf = new GurlsOptionsList("perf");
        opt.addOpt("perf", perf);

    //    if (pred_opt->getDataID() != typeid(T))
    //        return;

        gMat2D<T> tmp_pred(pred->cols(), pred->rows());

        const int pred_size = pred->getSize();
        //T* tmp_pred = new T[pred_size];

        Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));

        T* ap = new T[tot*t];


        //	for i = 1:tot
        for(int s = 0; s < tot; ++s)
        {
            //	guesses(i) = lmin*(q^i);
//            guesses[s] = lmin*std::pow(q,s+1);

            //		LL = L + (n*guesses(i));
            tmpvec = new T[xc];
            garbage.insert(tmpvec);
            set(tmpvec, n*guesses[s] , xc);
            axpy(xc, 1.0f, L, 1, tmpvec, 1);

            //		LL = LL.^(-1)
            setReciprocal(tmpvec, xc);
            //		LL = diag(LL);
            T* LL = diag(tmpvec, xc);
            garbage.insert(LL);

            delete[] tmpvec;
            garbage.erase(tmpvec);

            // WARNING: salvare tempo e memoria evitando di inizializzare
            // ad ogni passo

            //		num = y - LEFT*LL*RIGHT;
            tmp = new T[xc*t];
            garbage.insert(tmp);
            dot(LL, RIGHT, tmp, xc, xc, xc, t, xc, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

            tmp1 = new T[xr*t];
            garbage.insert(tmp1);
            dot(LEFT, tmp, tmp1, xr, xc, xc, t, xr, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

            delete[] tmp;
            garbage.erase(tmp);

            T* num = new T[xr*t];
            garbage.insert(num);
            // ?????
            copy(num, Y.getData(), xr*t);
            axpy(xr*t, -1.0f, tmp1, 1, num, 1);

//            copy(num, tmp1, xr*t);
//            axpy(xr*t, -1.0f, Y.getData(), 1, num, 1);

            delete[] tmp1;
            garbage.erase(tmp1);

            // den = zeros(n,1);
            set(den, 0.0f, n);

            //		for j = 1:n
            //			den(j) = 1-LEFT(j,:)*LL*right(:,j);
            //		end

            tmp = new T[xc];
            garbage.insert(tmp);
//			std::cout << std::endl << "got here!" << std::endl;

//			std::cout << LEFT.rows() << " - " << LEFT.cols() << std::endl;
//			std::cout << LL.rows() << " - " << LL.cols() << std::endl;
//			std::cout << tmp->rows() << " - " << tmp->cols() << std::endl;

            T* row = new T[xc];
            garbage.insert(row);

            for (int j = 0; j < n; ++j)
            {
                dot(LL, right +(xc*j), tmp, xc, xc, xc, 1, xc, 1, CblasNoTrans, CblasNoTrans, CblasColMajor);

                //extract j-th row from LEFT
                copy(row, LEFT + j, xc, 1, xr);

                den[j] =  ((T) 1.0) - dot (xc, row, 1, tmp, 1);
            }

            delete[] row;
            garbage.erase(row);
            delete[] tmp;
            garbage.erase(tmp);
            delete[] LL;
            garbage.erase(LL);



    //        opt.pred = zeros(n,T);
            set(tmp_pred.getData(), (T)0.0, pred_size);

            T* num_div_den = new T[n];

    //        for t = 1:T
            for(int j = 0; j< t; ++j)
            {
                rdivide(num + (n*j), den, num_div_den, n);

    //            opt.pred(:,t) = y(:,t) - (num(:,t)./den);
                minus(Y.getData() + (n*j), num_div_den, tmp_pred.getData()+(n*j), n);
            }

            delete [] num_div_den;
            tmp_pred.transpose(*pred);

    //        opt.perf = opt.hoperf([],y,opt);
            const gMat2D<T> dummy;
            perfClass->execute(dummy, Y_OMR, opt);

            gMat2D<T> *forho_vec = &(OptMatrix<gMat2D<T> >::dynacast(perf->getOpt("forho")))->getValue();

    //        for t = 1:T
            for(int j = 0; j<t; ++j)
            {
    //            ap(i,t) = opt.perf.forho(t);
                ap[s +(tot*j)] = forho_vec->getData()[j];
            }

            delete[] num;
            garbage.erase(num);
        }

        delete[] L;
        garbage.erase(L);
        delete[] LEFT;
        garbage.erase(LEFT);
        delete[] RIGHT;
        garbage.erase(RIGHT);
        delete[] right;
        garbage.erase(right);
        delete[] den;
        garbage.erase(den);
//        delete[] Le;
//        garbage.erase(Le);


        //[dummy,idx] = max(ap,[],1);
        const unsigned int* idx = indicesOfMax(ap, tot, t, 1);

        //vout.lambdas = 	guesses(idx);
        T* lambdas = copyLocations(idx, guesses, t, tot);

        delete[] idx;

        OptNumberList* LAMBDA = new OptNumberList();
        for (T* l_it = lambdas, *l_end = lambdas+t; l_it != l_end; ++l_it)
            LAMBDA->add(static_cast<double>(*l_it));

        delete[] lambdas;

        opt.addOpt("lambdas", LAMBDA);

        //vout.looe{1} = 	ap;

        tmp = transpose_rm(ap, t, tot);

        delete[] ap;

        gMat2D<T>* looe_mat = new gMat2D<T>(tmp, tot, t, true);

        delete[] tmp;

        opt.addOpt("looe", new OptMatrix<gMat2D<T> >(*looe_mat));


        //vout.guesses = 	guesses;
        gVec<T> guessesVector(guesses, tot, false);
        opt.addOpt("guesses", new OptMatrix<gMat2D<T> >(guessesVector.asMatrix()));

        delete[] guesses;

        opt.removeOpt("pred");
        if(pred_old != NULL)
            opt.addOpt("pred", pred_old);

    }
    catch( gException& e)
    {

        for(typename std::set<T*>::iterator it = garbage.begin(); it != garbage.end(); ++it)
            delete[] (*it);

        throw e;
    }
}

}

#endif // _GURLS_LOOCVPRIMAL_H_
