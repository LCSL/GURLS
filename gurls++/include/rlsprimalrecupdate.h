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


#ifndef _GURLS_RLSPRIMALRECUPDATE_H_
#define _GURLS_RLSPRIMALRECUPDATE_H_

#include "optimization.h"

#include "optmatrix.h"
#include "optfunction.h"

#include "utils.h"

namespace gurls
{

/**
 * \ingroup Optimization
 * \brief RLSPrimalRecUpdate is the sub-class of Optimizer that implements RLS with the primal formulation
 */
template <typename T>
class RLSPrimalRecUpdate: public Optimizer<T>
{
public:
    /**
     * Computes a classifier for the primal formulation of RLS.
     * The regularization parameter is set to the one found in the field paramsel of opt.
     * In case of multiclass problems, the regularizers need to be combined with the function specified inthe field singlelambda of opt
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - singlelambda (default)
     *  - paramsel (settable with the class ParamSelection and its subclasses)
     *
     * \return adds to opt the field optimizer which is a list containing the following fields:
     *  - W = matrix of coefficient vectors of rls estimator for each class
     *  - C = empty matrix
     *  - X = empty matrix
     *  - Cinv = inverse of the regularized kernel matrix in the primal space
     *  - XtX = the kernel matrix in the primal space
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* RLSPrimalRecUpdate<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //	[n,d] = size(X);

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();

    const unsigned long t = Y.cols();

    //  W = opt.rls.W;
    const gMat2D<T>& prev_W = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.W");
    gMat2D<T>* W = new gMat2D<T>(prev_W);

    //  Cinv = opt.rls.Cinv;
    const gMat2D<T>& prev_Cinv = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.Cinv");
    gMat2D<T>* Cinv = new gMat2D<T>(prev_Cinv);

    T* WData = W->getData();
    const unsigned long wn = W->rows();
    const unsigned long wd = W->cols();


    T* CinvData = Cinv->getData();
    const unsigned long cn = Cinv->rows();
    const unsigned long cd = Cinv->cols();

    T* Cx = new T[cn];
    T* x = new T[d];
    T* y = new T[t];
    T xCx;
    T* CxCxt = new T[cn*cn];
    T* xW = new T[wd];
    T* Cxy = new T[cn*t];

    for(unsigned long i=0; i<n; ++i)
    {
        getRow(X.getData(), n, d, i, x);
        getRow(Y.getData(), n, t, i, y);

        //  Cx = Cinv*X(i,:)';
        gemv(CblasNoTrans, cn, cd, (T)1.0, CinvData, cn, x, 1, (T)0.0, Cx, 1);

        //  xCx = X(i,:)*Cx;
        gemv(CblasNoTrans, 1, d, (T)1.0, x, 1, Cx, 1, (T)0.0, &xCx, 1);

        //  Cinv = Cinv - Cx*Cx'./(1+xCx);
        dot(Cx, Cx, CxCxt, cn, 1, cn, 1, cn, cn, CblasNoTrans, CblasTrans, CblasColMajor);
        axpy(cn*cn, (T)(-1.0/(xCx+1)), CxCxt, 1, CinvData, 1);


        //  W = W +(Cx*(y(i,:)-X(i,:)*W))./(1+xCx);
        dot(x, WData, xW, 1, d, wn, wd, 1, wd, CblasNoTrans, CblasNoTrans, CblasColMajor);
        axpy(t, (T)-1.0, xW, 1, y, 1);
        dot(Cx, y, Cxy, cn, 1, 1, t, cn, t, CblasNoTrans, CblasNoTrans, CblasColMajor);
        axpy(cn*t, (T)(1.0/(xCx+1)), Cxy, 1, WData, 1);
    }

    delete[] Cx;
    delete[] x;
    delete[] y;
    delete[] CxCxt;
    delete[] xW;
    delete[] Cxy;


    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

    //  rls.W = W;
    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    //  rls.C = [];
    gMat2D<T>* emptyC = new gMat2D<T>();
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

    //	cfr.X = [];
    gMat2D<T>* emptyX = new gMat2D<T>();
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));

    //  rls.Cinv = Cinv;
    optimizer->addOpt("Cinv", new OptMatrix<gMat2D<T> >(*Cinv));

    return optimizer;
}


}
#endif // _GURLS_RLSPRIMALRECUPDATE_H_
