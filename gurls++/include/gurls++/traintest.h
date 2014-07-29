/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
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

#ifndef GURLS_TRAINTEST_H
#define GURLS_TRAINTEST_H

#include <iostream>
#include <string>

#include "gurls++/gvec.h"
#include "gurls++/gmat2d.h"
#include "gurls++/optlist.h"
#include "gurls++/wrapper.h"
#include "gurls++/kernelrlswrapper.h"

namespace gurls{

//"easy train" fucntion
//requires in input:
//		- X, a T(float or double) buffer, intended as a matrix n x d
//		- Y, a T(float or double) buffer, intended as a matrix n x t
//		- n, d and t, geometry parameters described above
//		- algorithm, an optional string containing the selected algorithm (only "krls" available now)
//		- kernel, an optional string containing the selected kernel type ("gaussian" and "linear" available now)
//		- problem, an optional string containing the selected problem type ("classification" and "regression" available)
//		- savefile, an optional string containing the path in which to save the model
//gives in output
//		- a GurlsOptionList pointer containing the trained model, ready to be used with "test" function
template <typename T>
GurlsOptionsList* train(T* X, T* y, unsigned long n, unsigned long d, unsigned long t, 
		   std::string algorithm="krls", std::string kerneltype="gaussian", std::string problem="", std::string savefile="");

//"easy test" fucntion
//requires in input:
//		- model, previously saved opt, already trained on data using "train" function
//		- X, a T(float or double) buffer, intended as a matrix n x d
//		- predbuff, a T(float or double) buffer, intended as a matrix n x t, where prediction will be saved
//		- Y, an optional T(float or double) buffer, intended as a matrix n x t
//		- perfbuff, an optional T(float or double) buffer, intended as a t length vector, where performance will be saved
//		- n, d and t, geometry parameters described above
//		- perfstring, a string containing the type of perf to be calculated, "rmse", "macroavg" and "auto" are the accepted values,
//					  leave blank to skip performance calculation (perfbuff will remain untouched)
//gives in output
//		- EXIT_SUCCESS or EXIT_FAILURE
template <typename T>
int test(gurls::GurlsOptionsList& model, T* X, T* Y, T* predbuff, T* perfbuff, unsigned long n, unsigned long d, unsigned long t, std::string perfstring="");

//"easy test" function, with load from file
//requires in input:
//		- loadfile, path to previously saved optionlist in .bin format, already trained on data using "train" function
//		- X, a T(float or double) buffer, intended as a matrix n x d
//		- predbuff, a T(float or double) buffer, intended as a matrix n x t, where prediction will be saved
//		- Y, an optional T(float or double) buffer, intended as a matrix n x t
//		- perfbuff, an optional T(float or double) buffer, intended as a t length vector, where performance will be saved
//		- n, d and t, geometry parameters described above
//		- perfstring, a string containing the type of perf to be calculated, "rmse", "macroavg" and "auto" are the accepted values,
//					  leave blank to skip performance calculation (perfbuff will remain untouched)
//gives in output
//		- EXIT_SUCCESS or EXIT_FAILURE
template <typename T>
int test(std::string loadfile, T* X, T* Y, T* predbuff, T* perfbuff, unsigned long n, unsigned long d, unsigned long t, std::string perfstring="");
}

#include"traintest.hpp"

#endif
