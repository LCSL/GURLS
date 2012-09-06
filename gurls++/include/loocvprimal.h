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
#include "rmse.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamSelLoocvPrimal is the sub-class of ParamSelection that implements LOO cross-validation with the primal formulation
 */

template <typename T>
class ParamSelLoocvPrimal: public ParamSelection<T>
{

public:
    /**
     * Performs parameter selection when the primal formulation of RLS is used.
     * The leave-one-out approach is used.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following default fields:
     *  - nlambda (default)
     *  - smallnumber
     * \return adds the following fields to opt:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - guesses = array of guesses for the regularization parameter lambda
     *  - acc = matrix of validation accuracies for each lambda guess and for each class
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList *ParamSelLoocvPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    typename std::set<T*> garbage;

    try
    {
        //[n,T]  = size(y);
        const unsigned long n = Y.rows();
        const unsigned long t = Y.cols();


        //	K = X'*X;
        const unsigned long xc = X.cols();
        const unsigned long xr = X.rows();

        if(xr != n)
            throw gException(Exception_Inconsistent_Size);

        T* K = new T[xc*xc];
        garbage.insert(K);
        dot(X.getData(), X.getData(), K, xr, xc, xr, xc, xc, xc, CblasTrans, CblasNoTrans, CblasColMajor);


        // tot = opt.nlambda;
        int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));

        //	[Q,L] = eig(K);
        //	L = diag(L);

//        T* Q, *L;

//        eig(K, Q, L, xc, xc);
//        garbage.insert(Q);
//        garbage.insert(L);

        T* Q = K;
        T* L = new T[xc];
        garbage.insert(L);

        eig_sm(Q, L, xc);


        T* filtered = L;
        T* guesses = lambdaguesses(filtered, xc, std::min(xc,xr), xr, tot, (T)(opt.getOptAsNumber("smallnumber")));
        garbage.insert(guesses);
//        set(guesses, 0.f, tot);

        T* LOOSQE = new T[tot*t];
        garbage.insert(LOOSQE);
        set(LOOSQE, (T)0.0, tot*t);

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


        GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");

        gMat2D<T>* pred = new gMat2D<T>(n, t);
        OptMatrix<gMat2D<T> >* pred_opt = new OptMatrix<gMat2D<T> >(*pred);
        nestedOpt->addOpt("pred", pred_opt);

//        const int pred_size = pred->getSize();

        Performance<T>* perfClass = Performance<T>::factory(opt.getOptAsString("hoperf"));

        gMat2D<T>* perf = new gMat2D<T>(tot, t);
        T* ap = perf->getData();

        //	for i = 1:tot
        for(int s = 0; s < tot; ++s)
        {

            //		LL = L + (n*guesses(i));
            tmpvec = new T[xc];
            garbage.insert(tmpvec);
            set(tmpvec, n*guesses[s] , xc);
            axpy(xc, (T)1.0, L, 1, tmpvec, 1);

            //		LL = LL.^(-1)
            setReciprocal(tmpvec, xc);
            //		LL = diag(LL);
            T* LL = diag(tmpvec, xc);
            garbage.insert(LL);

            delete[] tmpvec;
            garbage.erase(tmpvec);

            // TODO WARNING: salvare tempo e memoria evitando di inizializzare
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

            copy(num, Y.getData(), xr*t);
            axpy(xr*t, (T)-1.0, tmp1, 1, num, 1);


            delete[] tmp1;
            garbage.erase(tmp1);

            // den = zeros(n,1);
            set(den, (T)0.0, n);

            //		for j = 1:n
            //			den(j) = 1-LEFT(j,:)*LL*right(:,j);
            //		end

            tmp = new T[xc];
            garbage.insert(tmp);

            T* row = new T[xc];
            garbage.insert(row);

            for (unsigned long j = 0; j < n; ++j)
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
//            set(pred->getData(), (T)0.0, pred_size);

            T* num_div_den = new T[n];

    //        for t = 1:T
            for(unsigned long j = 0; j< t; ++j)
            {
                rdivide(num + (n*j), den, num_div_den, n);

    //            opt.pred(:,t) = y(:,t) - (num(:,t)./den);
                copy(pred->getData()+(n*j), Y.getData() + (n*j), n);
                axpy(n, (T)-1.0, num_div_den, 1, pred->getData()+(n*j), 1);
            }

            delete [] num_div_den;

    //        opt.perf = opt.hoperf([],y,opt);
            const gMat2D<T> dummy;
            GurlsOptionsList* perf = perfClass->execute(dummy, Y, *nestedOpt);

            gMat2D<T> &forho_vec = perf->getOptValue<OptMatrix<gMat2D<T> > >("forho");

    //        for t = 1:T
            copy(ap+s, forho_vec.getData(), forho_vec.getSize(), tot, 1);

            delete perf;
            delete[] num;
            garbage.erase(num);
        }

        delete nestedOpt;
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
        unsigned long* idx = new unsigned long[t];
        T* work = NULL;
        indicesOfMax(ap, tot, t, idx, work, 1);

        //vout.lambdas = 	guesses(idx);
        T* lambdas = copyLocations(idx, guesses, t, tot);

        delete[] idx;


        gMat2D<T> *LAMBDA = new gMat2D<T>(1, t);
        copy(LAMBDA->getData(), lambdas, t);

        delete[] lambdas;


        GurlsOptionsList* paramsel;

        if(opt.hasOpt("paramsel"))
        {
            GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
            tmp_opt->copyOpt<T>("paramsel", opt);

            paramsel = GurlsOptionsList::dynacast(tmp_opt->getOpt("paramsel"));
            tmp_opt->removeOpt("paramsel", false);
            delete tmp_opt;

            paramsel->removeOpt("guesses");
            paramsel->removeOpt("perf");
            paramsel->removeOpt("lambdas");
        }
        else
            paramsel = new GurlsOptionsList("paramsel");



        paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*LAMBDA));
        paramsel->addOpt("perf", new OptMatrix<gMat2D<T> >(*perf));

        //vout.guesses = 	guesses;
        gMat2D<T> *guesses_mat = new gMat2D<T>(guesses, 1, tot, true);
        paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));

        delete[] guesses;

        return paramsel;
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
