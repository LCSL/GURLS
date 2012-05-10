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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cassert>

#include "gmat2d.h"
#include "gmath.h"
#include "options.h"
#include "optlist.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/base_object.hpp>

using namespace gurls;
using namespace std;

typedef float T;
const int Ncols = 12;
const int Nrows = 6;
const int Nsquared = 10;
const T tolerance = 1e-4;

#ifdef  USE_BINARY_ARCHIVES
typedef boost::archive::binary_iarchive iarchive;
typedef boost::archive::binary_oarchive oarchive;
#else
typedef boost::archive::text_iarchive iarchive;
typedef boost::archive::text_oarchive oarchive;
#endif

int main(int argc, char *argv[])
{

	cout.precision(4);
	cout.width(6);
	cout << fixed << right;
	try {
		std::cout << "Opening input and output files...";
		std::ofstream ofs("foo.txt", ios_base::binary);
		oarchive oa(ofs);
		std::ifstream ifs("foo.txt", ios_base::binary);
		std::cout << " done." << std::endl;

		std::cout << "Creating matrices and vectors...";
		gMat2D<T> A(Nrows,Ncols), B(Nsquared, Nsquared);
		gVec<T> V(Nsquared), U(Nsquared);

		B.randomize();
		A = 2;
		V = 2;
		U = 0;
		assert(A.sum()!=B.sum());
		assert(U.sum()!=V.sum());
		std::cout << " done." << std::endl;

		std::cout << "Saving data...";
		oa << B;
		oa << V;
		ofs.close();
		std::cout << " done." << std::endl;

		std::cout << "Loading data...";
		iarchive ia(ifs);
		ia >> A;
		ia >> U;
		ifs.close();
		std::cout << " done." << std::endl;

		std::cout << "Checking data integrity...";
		assert(A.sum()==B.sum());
		assert(U.sum()==V.sum());
		std::cout << " done." << std::endl;

		std::cout << std::endl;
		std::cout << "    *°*°*°*°*°*°*°*°*    " << std::endl;
		std::cout << std::endl;

		bool testString = false;
		bool testNumber = false;
		bool testMatrix = false;
		bool testTasks = false;
		bool testOptList = true;

		if (testString)  {
			OptString s("ciao");
			std::ofstream oparstream("par.txt");
			oarchive oparar(oparstream);
			oparar << s;
			oparstream.close();
			std::ifstream iparstream("par.txt");
			iarchive iparar(iparstream);
			OptString s1("pippo");
			std::cout << s << " - " << s1 << endl;
			iparar >> s1;
			iparstream.close();
			std::cout << s << " - " << s1 << endl;
		} else if (testNumber) {
			std::vector<string> sv;
			sv.push_back("good ");
			sv.push_back("night ");
			sv.push_back("and ");
			sv.push_back("good ");
			sv.push_back("luck!");
			OptStringList s(sv);
			std::ofstream oparstream("par.txt");
			oarchive oparar(oparstream);
			oparar << s;
			oparstream.close();
			std::ifstream iparstream("par.txt");
			iarchive iparar(iparstream);
			std::vector<string> sv1;
			sv1.push_back("good ");
			sv1.push_back("morning ");
			sv1.push_back("dude!");
			OptStringList s1(sv1);
			std::cout << s << " - " << s1 << endl;
			iparar >> s1;
			iparstream.close();
			std::cout << s << " - " << s1 << endl;
		} else if (testMatrix) {
			gMat2D<T> A(Nrows,Ncols);
			A = 1.5;
			OptMatrix<gMat2D<T> > s(A);
			std::ofstream oparstream("par.txt");
			oarchive oparar(oparstream);
			oparar << s;
			oparstream.close();
			std::ifstream iparstream("par.txt");
			iarchive iparar(iparstream);
			gMat2D<T> B(Ncols, Nrows);
			B.randomize();
			OptMatrix<gMat2D<T> > s1(B);
			std::cout << s << " - " << s1 << endl;
			iparar >> s1;
			iparstream.close();
			std::cout << s << " - " << s1 << endl;
		} else if (testTasks){

			OptTaskSequence s;
			s.addTask("paramsel:loocvprimal");
			s.addTask("optimizer:rlsprimal");
			s.addTask("pred:rlsprimal");
			s.addTask("perf:precrec");
			s.addTask("perf:macroavg");
			std::ofstream oparstream("par.txt");
			oarchive oparar(oparstream);
			oparar << s;
			oparstream.close();
			std::ifstream iparstream("par.txt");
			iarchive iparar(iparstream);
			OptTaskSequence s1;
			s1.addTask("paramsel:loocvprimal");
			std::cout << s << " - " << s1 << endl;
			iparar >> s1;
			iparstream.close();
			std::cout << s << " - " << s1 << endl;
		} else if (testOptList) {

			GurlsOptionsList *s = new GurlsOptionsList("prova", true);
			std::cout << *s << std::endl;

			std::ofstream oparstream("par.txt");
			oarchive oparar(oparstream);
			oparar << *s;
			oparstream.close();
			std::ifstream iparstream("par.txt");
			iarchive iparar(iparstream);
			GurlsOptionsList *s1 = new GurlsOptionsList("prova1", false);
			std::cout << *s1 << std::endl;
			iparar >> *s1;
			iparstream.close();
			std::cout << *s1 << std::endl;
		}


		std::cout << std::endl;
		std::cout << "    *°*°*°*°*°*°*°*°*    " << std::endl;
		std::cout << std::endl;

	}
	catch (gException& e){
		cout << e.getMessage() << endl;
		return -1;
	}

	return 0;

}
