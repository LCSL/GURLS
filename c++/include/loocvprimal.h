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

#include "options.h"
#include "optlist.h"
#include "gmat2d.h"
#include "gvec.h"
#include "gmath.h"

#include "paramsel.h"

namespace gurls {

template <typename T>
class LoocvPrimal: public ParamSelection<T>{

public:
	void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void LoocvPrimal<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt){

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
			guesses[s] = lmin*std::pow(q,s);

			//		LL = L + (n*guesses(i));
			//		LL = LL.^(-1);
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

}

#endif // _GURLS_LOOCVPRIMAL_H_
