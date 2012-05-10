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


#ifndef _GURLS_PRED_H_
#define _GURLS_PRED_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <exception>
#include <stdexcept>

#include "gmath.h"
#include "options.h"
#include "optlist.h"

using namespace std;

namespace gurls {

template <typename T>
class PredPrimal;

template <typename T>
class PredDual;

template <typename T>
class Prediction
{
public:
    virtual void execute( const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt) = 0;

    class BadPredictionCreation : public std::logic_error {
    public:
        BadPredictionCreation(std::string type)
            : logic_error("Cannot create type " + type) {}
    };

    static Prediction<T>*
    factory(const std::string& id) throw(BadPredictionCreation) {
        if(id == "rlsprimal"){
            return new PredPrimal<T>;
        }	else if(id == "rlsdual"){
            return new PredDual<T>;
        } else
            throw BadPredictionCreation(id);
    }

};

//template <typename Matrix>
//class PredPrimal: public Prediction<Matrix> {

//public:
//	Matrix& execute(const Matrix &X, const Matrix &Y, GurlsOptionsList &opt);
//};

//template <typename Matrix>
//class PredDual: public Prediction<Matrix> {

//public:
//	Matrix& execute(const Matrix &X, const Matrix &Y, GurlsOptionsList& opt);
//};

//template <typename Matrix>
//Matrix& PredPrimal<Matrix>::execute(const Matrix &X, const Matrix &Y,
//									GurlsOptionsList &opt){
//	if (opt.hasOpt("W")){
//		GurlsOption *g = opt.getOpt("W");
//		//Matrix& W = OptMatrixTMP<Matrix>::dynacast(g)->getValue();
//		Matrix& W = OptMatrix<Matrix>::dynacast(g)->getValue();
//		Matrix* Z = new Matrix(X.rows(), W.cols());
//		*Z = 0;
//		dot(X, W, *Z);
//		return *Z;
//	}else {
//		throw gException(gurls::Exception_Required_Parameter_Missing);
//	}
//}

//template <typename Matrix>
//Matrix& PredDual<Matrix>::execute(const Matrix& X, const Matrix &Y,
//								  GurlsOptionsList& opt){

//	try {

//		string type = OptString::dynacast( opt.getOpt("kernel.type") )->getValue();
//		if (type == "linear"){
//			PredPrimal<Matrix> pred;
//			return pred.execute(X, Y, opt);

//		}else {
//			if (opt.hasOpt("finalkernel.K") && opt.hasOpt("rls.C")){
//				GurlsOption *g = opt.getOpt("finalkernel.K");
//				//Matrix& K = OptMatrixTMP<Matrix>::dynacast(g)->getValue();
//				Matrix& K = OptMatrix<Matrix>::dynacast(g)->getValue();
//				g = opt.getOpt("rls.C");
//				//Matrix& C = OptMatrixTMP<Matrix>::dynacast(g)->getValue();
//				Matrix& C = OptMatrix<Matrix>::dynacast(g)->getValue();
//				Matrix* Z = new Matrix(K->rows(), C->cols());
//				*Z = 0;
//				dot(*K, *C, *Z);
//				return *Z;
//			}else {
//				throw gException(gurls::Exception_Required_Parameter_Missing);
//			}
//		}
//	}catch (gException& ex){
//		throw ex;
//	}

//}

}
#endif // _GURLS_PRED_H_
