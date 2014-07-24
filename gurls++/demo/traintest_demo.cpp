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

		if (argc<2){
			std::cout<<"\tUsage:\t traintest_demo \"data folder\""<<std::endl;
			return EXIT_FAILURE;}
		gMat2D<T> Xtr, ytr, Xte, yte;
		try
		{
			std::cout<<"Reading Xtr"<<std::endl;
        Xtr.readCSV(std::string(argv[1]) + "/Xtr.txt");
			std::cout<<"Reading ytr"<<std::endl;
        ytr.readCSV(std::string(argv[1]) + "/ytr.txt");
			std::cout<<"Reading Xte"<<std::endl;
        Xte.readCSV(std::string(argv[1]) + "/Xte.txt");
			std::cout<<"Reading yte"<<std::endl;
        yte.readCSV(std::string(argv[1]) + "/yte.txt");
		}
		catch(gurls::gException &e)
		{
			std::cout<<e.getMessage();
			return EXIT_FAILURE;
		}

		T* perfBuffer = new T[yte.cols()];
		T* predBuffer = new T[yte.cols()*yte.rows()];

		std::cout<<"Testing with (kernel), keeping opt in memory"<<std::endl;	
		try
		{
		gurls::GurlsOptionsList opt = train(Xtr.getData(), ytr.getData(), Xtr.rows(), Xtr.cols(), ytr.cols(), "krls", "gaussian");
		test(opt, Xte.getData(), yte.getData(), predBuffer, perfBuffer, Xte.rows(), Xte.cols(), yte.cols(), "auto");
		}
		catch(gException &e)
		{
		std::cout<<e.what();
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
		
		std::cout<<"Testing with (linear), writing opt to file"<<std::endl;
		try
		{
		train(Xtr.getData(), ytr.getData(), Xtr.rows(), Xtr.cols(), ytr.cols(), "krls", "linear","","tempfile");
		test("tempfile", Xte.getData(), yte.getData(), predBuffer, perfBuffer, Xte.rows(), Xte.cols(), yte.cols(), "auto");
		std::remove("tempfile");
		}
		catch(gException &e)
		{
		std::cout<<e.what();
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
