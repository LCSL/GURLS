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

#include <iostream>
#include <string>
#include "gurls++/gurls.h"
#include "gurls++/traintest.h"

using namespace gurls;

int main(int argc, char* argv[])
{
		typedef double T;

		if (argc!=2){
			std::cout<<"\tUsage:\t traintest_demo \"data folder\""<<std::endl;
			return EXIT_FAILURE;}
		//gurls Matrix objects
		gMat2D<T> Xtr, ytr, Xte, yte;

		try
		{
		//reading CSV with standard names from selected path
			std::cout<<"Reading Xtr"<<std::endl;
        Xtr.readCSV(std::string(argv[1]) + "/Xtr.csv");
			std::cout<<"Reading ytr"<<std::endl;
        ytr.readCSV(std::string(argv[1]) + "/Ytr.csv");
			std::cout<<"Reading Xte"<<std::endl;
        Xte.readCSV(std::string(argv[1]) + "/Xts.csv");
			std::cout<<"Reading yte"<<std::endl;
        yte.readCSV(std::string(argv[1]) + "/Yts.csv");
		}

		catch(gurls::gException &e)
		{
			std::cout<<e.getMessage();
			return EXIT_FAILURE;
		}
		
		//allocating the buffers where to write prediction and performance information
		//perfBuffer must be greater or equal than the number of columns of Y
		T* perfBuffer = new T[yte.cols()];
		//predBuffer must be greater or equal to total size of Yts
		T* predBuffer = new T[yte.cols()*yte.rows()];

		std::cout<<"Running train and test functions..."<<std::endl;	
		try
		{
		unsigned long n= Xtr.rows();
		unsigned long d= Xtr.cols();
		unsigned long t= ytr.cols();
		unsigned long nTest= Xte.rows();

		//train function returns a GurlsOptionsList pointer, containing the trained model parameters, "linear" is selected as kernel value in this demo
		gurls::GurlsOptionsList* model = train(Xtr.getData(), ytr.getData(), n, d, t, "krls", "linear");
		//test writes directly prediction and performance in predBuffer and perfBuffer, "auto" is selected as performance type in this demo
		test(*model, Xte.getData(), yte.getData(), predBuffer, perfBuffer, nTest, d, t, "auto");
		
		delete model;
		}
		catch(gException &e)
		{
			std::cout<<e.getMessage();
			delete [] perfBuffer;
			delete [] predBuffer;
			return EXIT_FAILURE;
		}
		catch(std::exception &e)
		{
			std::cout<<e.what();
			delete [] perfBuffer;
			delete [] predBuffer;
			return EXIT_FAILURE;
		}

		std::cout<<"Performance:"<<std::endl;
		for(unsigned int i=0; i<yte.cols(); ++i)
			std::cout<<perfBuffer[i]<<" ";
		std::cout<<std::endl;
		
		delete [] perfBuffer;
		delete [] predBuffer;
}
