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


#ifndef _GURLS_RLSPRIMAL_H_
#define _GURLS_RLSPRIMAL_H_

#include "optimization.h"

namespace gurls {

template <typename T>
class RLSPrimal: public Optimizer<T>{

public:
	void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};


template <typename T>
void RLSPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt){

	//	lambda = opt.singlelambda(opt.paramsel.lambdas);
	std::vector<double> ll = OptNumberList::dynacast(opt.getOpt("lambdas"))->getValue();
	double* lld  = new double[ll.size()];
	memcpy( lld, &ll[0], sizeof( double) * ll.size() );
	double lambda = (OptFunction::dynacast(opt.getOpt("singlelambda")))->getValue(lld, ll.size());

	//	fprintf('\tSolving primal RLS...\n');
	std::cout << "Solving primal RLS... ";
	//	[n,d] = size(X);
	T n = X.rows();
	T d = X.cols();

	//	===================================== Primal K
	//	K = X'*X;
	gMat2D<T>* K = new gMat2D<T>(X.cols(), X.cols());
	gMat2D<T>* Xt = new gMat2D<T>(X.cols(), X.rows());
	X.transpose(*Xt);
	dot(*Xt, X, *K);

	try{ // Try solving it with cholesky first.

		//		K = K + (n*lambda)*eye(d);
		*K += (n*static_cast<T>(lambda))*gMat2D<T>::eye(d);


		//		R = chol(K);
		gMat2D<T>* R = new gMat2D<T>(K->rows(), K->cols());
		cholesky(*K, *R);
		gMat2D<T>* Rt = new gMat2D<T>(R->cols(), R->rows());
		R->transpose(*Rt);
		//		cfr.dcomptime = toc;

		//		Xty = X'*y;
		gMat2D<T>* Xty = new gMat2D<T>(Xt->rows(), Y.cols());
		dot(*Xt, Y, *Xty);

		// if isfield(opt,'W0')
		if (opt.hasOpt("W0")){
			//			Xty = Xty + opt.W0;
			GurlsOption *g = opt.getOpt("W0");
			// check if the stored W0 has a compatible data type
			if (g->getDataID() == typeid(T)){
				//gMat2D<T> *W0 = &(OptMatrix<T>::dynacast(g))->getValue();
				gMat2D<T> *W0 = &(OptMatrix<gMat2D<T> >::dynacast(g))->getValue();
				*Xty += *W0;
			}
		}
		//		cfr.W = R\(R'\Xty);

		gMat2D<T>* R_inv = new gMat2D<T>(R->cols(), R->rows());
		pinv(*R, *R_inv);

		gMat2D<T>* Rt_inv = new gMat2D<T>(Rt->cols(), Rt->rows());
		pinv(*Rt, *Rt_inv);

		gMat2D<T>* Wtmp = new gMat2D<T>(Rt_inv->rows(), Xty->cols());
		dot(*Rt_inv, *Xty, *Wtmp);
		gMat2D<T>* W = new gMat2D<T>(R_inv->rows(), Wtmp->cols());
		dot(*R_inv, *Wtmp, *W);

		//opt.addOpt("W", new OptMatrixTMP<gMat2D<T> >(*W));
		opt.addOpt("W", new OptMatrix< gMat2D<T> > (*W));
		std::cout << "W = " << std::endl << *W << std::endl;
	} catch (gException& gex) {

		// AT PRESENT; THIS SECOND OPTION HAS NOT BEEN IMPEMENTED YET
		throw gex;
		//	catch
		//	%% If it fails, then solve it using svd and rls_eigen.
		//		[Q,L,V] = svd(K);
		//		Q = double(Q);
		//		L = double(diag(L));

		//		QtXtY = Q'*(X'*y);

		//		% regularization is done inside rls_eigen
		//		cfr.W = rls_eigen(Q, L, QtXtY, lambda, n);
		//	end

	}

	// WHAT SHALL WE DO WITH THE FOLLOWING TWO STATEMENTS?
	//	cfr.C = [];
	//	cfr.X = [];

}




}
#endif // _GURLS_RLSPRIMAL_H_
