 /*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
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


#ifndef _GURLS_RLSPRIMALRECINIT_H_
#define _GURLS_RLSPRIMALRECINIT_H_

#include "gurls++/optimization.h"

#include "gurls++/optmatrix.h"
#include "gurls++/optfunction.h"

#include "gurls++/utils.h"

namespace gurls
{

/**
 * \ingroup Optimization
 * \brief RLSPrimalRecInit is the sub-class of Optimizer that implements RLS with the primal formulation
 */
template <typename T>
class RLSPrimalRecInit: public Optimizer<T>
{

public:
	///
	/// Default constructor
	///
	RLSPrimalRecInit():Optimizer<T>("rlsprimalrecinit"){}
	
	///
	/// Clone method
	///
	TaskBase *clone()
	{
		return new RLSPrimalRecInit<T>();
	}

    /**
     * Computes a classifier for the primal formulation of RLS.
     * The regularization parameter is set to the one found in the field paramsel of opt.
	 * The variables necessary for further recursive update are stored in the output structure
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
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* RLSPrimalRecInit<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //	lambda = opt.singlelambda(opt.paramsel.lambdas);
    const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
    T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());


    //	[n,d] = size(X);
    const unsigned long n = opt.hasOpt("nTot")? static_cast<unsigned long>(opt.getOptAsNumber("nTot")) : X.rows();
    unsigned long d;
    unsigned long t;


    //	XtX = X'*X;
    T* XtX;
    if(!opt.hasOpt("kernel.XtX"))
    {
        d = X.cols();
        XtX = new T[d*d];
        dot(X.getData(), X.getData(), XtX, n, d, n, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);
    }
    else
    {
        const gMat2D<T>& XtX_mat = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.XtX");
        d = XtX_mat.cols();
        XtX = new T[d*d];
        copy(XtX, XtX_mat.getData(), d*d);
    }

    //	Xty = X'*y;
    T* Xty;
    if(!opt.hasOpt("kernel.Xty"))
    {
        t = Y.cols();
        Xty = new T[d*t];
        dot(X.getData(), Y.getData(), Xty, n, d, n, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
    }
    else
    {
        const gMat2D<T>& Xty_mat = opt.getOptValue<OptMatrix<gMat2D<T> > >("kernel.Xty");
        t = Xty_mat.cols();
        Xty = new T[d*t];
        copy(Xty, Xty_mat.getData(), d*t);
    }


    //  Cinv = pinv(XtX + (n*lambda)*eye(d));
    T coeff = n*lambda;
    axpy(d, (T)1.0, &coeff, 0, XtX, d+1);

    int cinv_rows, cinv_cols;
    T* cinv = pinv(XtX, d, d, cinv_rows, cinv_cols);

    delete[] XtX;

    gMat2D<T>*W = new gMat2D<T>(cinv_rows, t);
    dot(cinv, Xty, W->getData(), cinv_rows, d, d, t, cinv_rows, t, CblasNoTrans, CblasNoTrans, CblasColMajor);

    gMat2D<T> *Cinv = new gMat2D<T>(cinv_rows, cinv_cols);
    copy(Cinv->getData(), cinv, Cinv->getSize());

    delete[] cinv;
    delete[] Xty;

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

    //  cfr.W = Cinv*Xty;
    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    //  cfr.C = [];
    gMat2D<T>* emptyC = new gMat2D<T>();
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

    //	cfr.X = [];
    gMat2D<T>* emptyX = new gMat2D<T>();
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));

    //  cfr.Cinv = Cinv;
    optimizer->addOpt("Cinv", new OptMatrix<gMat2D<T> >(*Cinv));

    return optimizer;
}


}
#endif // _GURLS_RLSPRIMALRECINIT_H_
