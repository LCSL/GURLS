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


#ifndef _GURLS_BIGHOPRIMAL_H_
#define _GURLS_BIGHOPRIMAL_H_


#include "bigparamsel.h"
#include "optmatrix.h"
#include "bigmath.h"
#include "bigpred_primal.h"
#include "bigperf.h"


namespace gurls
{

/**
 * \ingroup ParameterSelection
 * \brief BigParamSelHoPrimal is the subclass of BigParamSelection that implements hold-out cross validation with the primal formulation of RLS
 */
template <typename T>
class BigParamSelHoPrimal: public BigParamSelection<T>
{

public:

    /**
     * Performs parameter selection when the primal formulation of RLS is used.
     * The hold-out approach is used.
     * The performance measure specified by opt.hoperf is maximized.
     * \param X input data bigarray
     * \param Y labels bigarray
     * \param opt options with the following:
     *  - nlambda (default)
     *  - hoperf (default)
     *  - smallnumber (default)
     *  - split (settable with the class Split and its subclasses)
     *
     * \return a GurlsOptionList with the following fields:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - guesses = array of guesses for the regularization parameter lambda
     *  - forho = matrix of validation accuracies for each lambda guess and for each class
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt);

};

template <typename T>
GurlsOptionsList *BigParamSelHoPrimal<T>::execute(const BigArray<T> &X, const BigArray<T> &Y, const GurlsOptionsList &opt)
{

    const GurlsOptionsList* split = opt.getOptAs<GurlsOptionsList>("split");
    const BigArray<T>& Xva = split->getOptValue<OptMatrix<BigArray<T> > >("Xva");
    const BigArray<T>& Yva = split->getOptValue<OptMatrix<BigArray<T> > >("Yva");

    BigArray<T>* XtX = matMult_AtB(X, X, opt.getOptAsString("files.XtX_filename"), opt.getOptAsNumber("memlimit"));
    BigArray<T>* Xty = matMult_AtB(X, Y, opt.getOptAsString("files.Xty_fileName"), opt.getOptAsNumber("memlimit"));
    BigArray<T>* XvatXva = matMult_AtB(Xva, Xva, opt.getOptAsString("files.XvatXva_filename"), opt.getOptAsNumber("memlimit"));
    BigArray<T>* Xvatyva = matMult_AtB(Xva, Yva, opt.getOptAsString("files.Xvatyva_fileName"), opt.getOptAsNumber("memlimit"));


//    K = XtX - XvatXva;
    gMat2D<T>* XtX_mat = new gMat2D<T>(XtX->rows(), XtX->cols());
    XtX->getMatrix(0,0, *XtX_mat);
    delete XtX;

    gMat2D<T>* XvatXva_mat = new gMat2D<T>(XvatXva->rows(), XvatXva->cols());
    XvatXva->getMatrix(0,0, *XvatXva_mat);
    delete XvatXva;

    axpy(XtX_mat->getSize(), (T)-1.0, XvatXva_mat->getData(), 1, XtX_mat->getData(), 1);
    delete XvatXva_mat;

//  Xty = Xty - Xvatyva;
    gMat2D<T>* Xty_mat = new gMat2D<T>(Xty->rows(), Xty->cols());
    Xty->getMatrix(0,0, *Xty_mat);
    delete Xty;

    gMat2D<T>* Xvatyva_mat = new gMat2D<T>(Xvatyva->rows(), Xvatyva->cols());
    Xvatyva->getMatrix(0,0, *Xvatyva_mat);
    delete Xvatyva;

    axpy(Xty_mat->getSize(), (T)-1.0, Xvatyva_mat->getData(), 1, Xty_mat->getData(), 1);
    delete Xvatyva_mat;


    const unsigned long n = X.rows() - Xva.rows();
    const unsigned long d = X.cols();
    const unsigned long t = Y.cols();


    int tot = static_cast<int>(std::ceil( opt.getOptAsNumber("nlambda")));


//	[Q,L] = eig(K);
    T* Q = XtX_mat->getData();
    T *L = new T[d];

    eig_sm(Q, L, d);


//	QtXtY = Q'*Xty;
    T* QtXtY = new T[d*t];
    dot(Q, Xty_mat->getData(), QtXtY, d, d, d, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
    delete Xty_mat;

//	guesses = paramsel_lambdaguesses(L, min(n,d), n, opt);
    T* guesses = lambdaguesses(L, d, std::min(d, n), n, tot, (T)(opt.getOptAsNumber("smallnumber")));

//	ap = zeros(tot,T);
    gMat2D<T>* ap_mat = new gMat2D<T>(tot, t);
    T* ap = ap_mat->getData();
    set(ap, (T)0.0, tot*t);


    GurlsOptionsList* nestedOpt = new GurlsOptionsList("nested");
    nestedOpt->copyOpt("nb_pred", opt);
    nestedOpt->copyOpt("files", opt);
    nestedOpt->copyOpt("tmpfile", opt);
    nestedOpt->copyOpt("memlimit", opt);

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");
    nestedOpt->addOpt("optimizer",optimizer);

    T* work = new T[d*(d+1)];

    gMat2D<T> *W = new gMat2D<T>(d, t);
    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    BigPredPrimal< T > primal;
    BigPerformance<T>* perfClass = BigPerformance<T>::factory(opt.getOptAsString("hoperf"));


    for(int i=0; i<tot; ++i)
    {
//        opt.rls.W = rls_eigen(Q,L,QtXtY,guesses(i),n);
        rls_eigen(Q, L, QtXtY, W->getData(), guesses[i], n, d, d, d, d, t, work);

//        opt.pred = bigpred_primal(Xva,yva,opt);
        GurlsOptionsList *ret_pred = primal.execute(Xva, Yva, *nestedOpt);

        nestedOpt->removeOpt("pred");
        nestedOpt->addOpt("pred", ret_pred);

//		opt.perf = opt.hoperf(Xva,yva,opt);
        GurlsOptionsList* ret_perf = perfClass->execute(Xva, Yva, *nestedOpt);

        gMat2D<T> &forho_vec = ret_perf->getOptValue<OptMatrix<gMat2D<T> > >("forho");

        copy(ap+i, forho_vec.getData(), t, tot, 1);

        delete ret_perf;
    }

    delete [] work;
    delete nestedOpt;

    delete perfClass;
    delete [] XtX_mat;
    delete [] L;
    delete [] QtXtY;


    GurlsOptionsList* paramsel;

    if(opt.hasOpt("paramsel"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt("paramsel", opt);

        paramsel = GurlsOptionsList::dynacast(tmp_opt->getOpt("paramsel"));
        tmp_opt->removeOpt("paramsel", false);
        delete tmp_opt;

        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("forho");
        paramsel->removeOpt("guesses");
    }
    else
        paramsel = new GurlsOptionsList("paramsel");



//    [dummy,idx] = max(ap,[],1);
    work = NULL;
    unsigned long* idx = new unsigned long[t];
    indicesOfMax(ap, tot, t, idx, work, 1);

//	vout.lambdas = guesses(idx);
    gMat2D<T>* lambdas = new gMat2D<T>(1, t);
    copyLocations(idx, guesses, t, tot, lambdas->getData());
    delete [] idx;

    paramsel->addOpt("lambdas", new OptMatrix<gMat2D<T> >(*lambdas));


//	vout.forho = ap;
    paramsel->addOpt("forho", new OptMatrix<gMat2D<T> >(*ap_mat));

//	vout.guesses = guesses;
    gMat2D<T>* guesses_mat = new gMat2D<T>(1, tot);
    copy(guesses_mat->getData(), guesses, tot);
    delete[] guesses;

    paramsel->addOpt("guesses", new OptMatrix<gMat2D<T> >(*guesses_mat));


    return paramsel;
}

}

#endif // _GURLS_BIGHOPRIMAL_H_
