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


#ifndef _GURLS_RLSPRIMAL_H_
#define _GURLS_RLSPRIMAL_H_

#include "gurls++/optimization.h"

#include "gurls++/optmatrix.h"
#include "gurls++/optfunction.h"

#include "gurls++/utils.h"

namespace gurls {

/**
 * \ingroup Optimization
 * \brief RLSPrimal is the sub-class of Optimizer that implements RLS with the primal formulation
 */
template <typename T>
class RLSPrimal: public Optimizer<T>{

public:	
	///
	/// Default constructor
	///
	RLSPrimal():Optimizer<T>("rlsprimal"){}
	
	///
	/// Clone method
	///
	TaskBase *clone()
	{
		return new RLSPrimal<T>();
	}

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
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* RLSPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList &opt)
{
    //	lambda = opt.singlelambda(opt.paramsel.lambdas);
    const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
    T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());

//    std::cout << "Solving primal RLS... " << std::endl;

    //	[n,d] = size(X);

    const long n = X.rows();
    const long d = X.cols();

    const long Yn = Y.rows();
    const long Yd = Y.cols();

    //	===================================== Primal K

    //	K = X'*X;
    T* K = new T[d*d];
    dot(X.getData(), X.getData(), K, n, d, n, d, d, d, CblasTrans, CblasNoTrans, CblasColMajor);

    //	Xty = X'*y;
    T* Xty = new T[d*Yd];
    dot(X.getData(), Y.getData(), Xty, n, d, Yn, Yd, d, Yd, CblasTrans, CblasNoTrans, CblasColMajor);


    gMat2D<T> *W = rls_primal_driver(K, Xty, n, d, Yd, lambda);

    delete[] K;
    delete[] Xty;

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    gMat2D<T>* emptyC = new gMat2D<T>();
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

    //	cfr.X = [];
    gMat2D<T>* emptyX = new gMat2D<T>();
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));

    return optimizer;
}


}
#endif // _GURLS_RLSPRIMAL_H_
