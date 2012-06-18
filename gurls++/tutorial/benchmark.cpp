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
#include <cstring>
#include <ctime>
#include <cassert>

#include "gvec.h"
#include "gmat2d.h"
#include "gmath.h"
#include "norm.h"

using namespace gurls;
using namespace std;

// Change this line to make a test with different data types
typedef float T;
const int CLOCKS_PER_MILLISEC = CLOCKS_PER_SEC/1000;
const int reportLength = 50;
const int TOT_COUNT = static_cast<int>(1E5);

void benchmarkGVEC(int N, int reps);
void benchmarkGMAT2D(int Nrows, int Ncols, int reps);
void benchmarkGMATH(int Nrows, int Ncols, int reps);
void benchmarkNORM(int Nrows, int Ncols, int reps);
void benchmarkKERNEL(int Nrows, int Ncols, int reps);

void printReportMsg(string operation, clock_t start, clock_t end, int reps){
	string filling(reportLength - operation.length(), ' ');
	cout << "[" << operation << filling << "]	required "
		 << (double(end - start) / CLOCKS_PER_MILLISEC / reps)
		 << " milliseconds on average." << endl;
}


int main(int argc, char *argv[])
{

	int N1 = static_cast<int>(1E4);
	int N2 = static_cast<int>(4E3);
	int reps = TOT_COUNT/N1;

	if (argc > 4){
		cout	<< "Error: too many arguments provided." << endl
				<< "Usage: " << argv[0] << " N1 N2 reps " << endl
				<< "where <N1> (optional) is either the length of the "
				<< "vectors or the number of rows of the matrices, <N2> "
				<< "(optional) is the number fo columns of the matrices, "
				<< "and <reps> (optional) is the number of repetitions of "
				<< "each test." << endl << endl;
		return -1;
	}

	if (argc > 1) {
		N1 = atoi( argv[1] );
		if (argc > 2){
			N2 = atoi( argv[2]);
		}
		if (N1 < TOT_COUNT){
			reps = static_cast<int>(ceil((double)( TOT_COUNT/N1)));
		} else {
			reps = 1;
		}
		if (argc > 3){
			reps = atoi( argv[3] );
		}
	}

	cout << endl;
	cout << "============================================= " << endl;
	cout << "========== STARTING GURLS BENCHMARK =========" << endl;
	cout << "============================================= " << endl;
	cout << "Nrows = " << N1 << endl << "Ncols = " << N2 << endl
		 << "Repetitions = "<< reps << endl << endl;

	cout << fixed << endl;
	cout.precision(6);

	cout << endl;
	cout << "============================================= " << endl;
	cout << "====== Testing methods of `gVec' class ====== " << endl;
	cout << "============================================= " << endl;
	benchmarkGVEC(N1, reps);

	cout << endl;
	cout << "============================================= " << endl;
	cout << "===== Testing methods of `gMat2D' class ===== " << endl;
	cout << "============================================= " << endl;
	benchmarkGMAT2D(N1, N2, reps);


	cout << endl;
	cout << "============================================= " << endl;
	cout << "===== Testing methods of `gMath' class ====== " << endl;
	cout << "============================================= " << endl;
	benchmarkGMATH(N1, N2, reps);

	cout << endl;
	cout << "============================================= " << endl;
	cout << "==== Testing methods defined in `norm.h' ==== " << endl;
	cout << "============================================= " << endl;
	benchmarkNORM(N1, N2, reps);




	cout << endl;
	cout << "============================================= " << endl;
	cout << "========== GURLS BENCHMARK COMPLETE =========" << endl ;
	cout << "============================================= " << endl << endl;

	return 0;
}

void benchmarkGVEC(int N, int reps){

	clock_t start,end;

	gVec<T> *A, *B, *C;
	try{
		// -------------------------------------------- CREATION AND ASSIGNMENT
		start = clock();
		for (int i = 0; i < reps; ++i) {
			gVec<T> a = gVec<T>::zeros(N);
		}
		end = clock();
		printReportMsg("vector creation and assignment", start, end, reps);

		// -------------------------------------------------- RANDOM ASSIGNMENT
		start = clock();
		for (int i = 0; i < reps; ++i) {
			gVec<T> a = gVec<T>::rand(N);
		}
		end = clock();
		printReportMsg("random vector creation", start, end, reps);

		// --------------------------------------------- ADDITION OF A CONSTANT
		A = new gVec<T>(N);
		start = clock();
		for (int i = 0; i < reps; ++i) {
			*A+10;
		}
		delete A;
		end = clock();
		printReportMsg("vector-constant addition", start, end, reps);

		// ----------------------------------------------------------- ADDITION
		A = new gVec<T>(N);
		B = new gVec<T>(N);
		start = clock();
		for (int i = 0; i < reps; ++i) {
			*A+*B;
		}
		end = clock();
		delete A;
		delete B;
		printReportMsg("vector-vector addition", start, end, reps);

		// ----------------------------------------------------- MULTIPLICATION
		A = new gVec<T>(N);
		B = new gVec<T>(N);
		start = clock();
		for (int i = 0; i < reps; ++i) {
			(*A) * (*B);
		}
		end = clock();
		delete A;
		delete B;
		printReportMsg("vector-vector multiplication", start, end, reps);

		// --------------------------------------------------- INPLACE DIVISION
		A = new gVec<T>(N);
		*A = TOT_COUNT;
		B = new gVec<T>(N);
		*B = 2;
		start = clock();
		for (int i = 0; i < reps; ++i) {
			(*A) /= (*B);
		}
		assert(A->at(0)*pow(B->at(0), reps) == TOT_COUNT);
		end = clock();
		delete A;
		delete B;
		printReportMsg("inplace vector-vector division", start, end, reps);

		// ---------------------------------------------- OUT-OF-PLACE DIVISION
		A = new gVec<T>(N);
		*A = 2;
		B = new gVec<T>(N);
		*B = 2;
		C = new gVec<T>(N);
		start = clock();
		for (int i = 0; i < reps; ++i) {
			*C = (*A) / (*B);
		}
		assert(C->at(0)==1);
		end = clock();
		assert( start != -1);
		assert( end != -1);
		delete A;
		delete B;
		delete C;
		printReportMsg("out-of-place vector-vector division", start, end, reps);


		// --------------------------------------------------- CHECK "CLOSE TO"
		A = new gVec<T>(N);
		B = new gVec<T>(N);
		*A = 2.0;
		*B = gVec<T>::rand(N)/static_cast<T>(1E11+2);
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A->closeTo(*B, (T)1E-10);
		}
		assert(A->closeTo(*B, (T)1E-10));
		end = clock();
		delete A;
		delete B;
		printReportMsg("test of method closeto()", start, end, reps);
	}
	catch (gException& e) {

		cout << e.getMessage() << endl;
	}

}


void benchmarkGMAT2D(int Nrows, int Ncols, int reps){


	clock_t start,end;

	gMat2D<T> A(Nrows,Ncols), B(Nrows,Ncols), C(Nrows,Ncols);
	try {
		// -------------------------------------------- CREATION AND ASSIGNMENT
		start = clock();
		for (int i = 0; i < reps; ++i) {
			gMat2D<T> a = gMat2D<T>::zeros(Nrows,Ncols);
		}
		end = clock();
		printReportMsg("matrix creation and assignment", start, end, reps);

		// -------------------------------------------------- RANDOM ASSIGNMENT
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A.randomize();
		}
		end = clock();
		printReportMsg("random matrix creation", start, end, reps);

		// --------------------------------------------- ADDITION OF A CONSTANT
		A = 0;
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A+10;
		}
		end = clock();
		printReportMsg("matrix-constant addition", start, end, reps);

		// ---------------------------------------------- OUT-OF-PLACE ADDITION
		A = 1.0;
		B = 1.0;
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A+B;
		}
		end = clock();
		printReportMsg("out-of-place matrix-matrix addition", start, end, reps);

		// -------------------------------------------------- IN-PLACE ADDITION
		A = 1.0;
		B = 1.0;
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A+=B;
		}
		end = clock();
		printReportMsg("inplace matrix-matrix addition", start, end, reps);

		// ---------------------------------------- OUT-OF-PLACE MULTIPLICATION
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A * B;
		}
		end = clock();
		printReportMsg("out-of-place matrix-matrix multiplication", start,
					   end, reps);

		// --------------------------------------------------- INPLACE DIVISION
		A = static_cast<T>(TOT_COUNT);
		B = 2.0;
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A /= B;
		}
		assert(A(0,0)*pow(B(0,0), reps) == TOT_COUNT);
		end = clock();
		printReportMsg("inplace matrix-matrix division", start, end, reps);

		// ---------------------------------------------- OUT-OF-PLACE DIVISION
		A = 2.0;
		B = 2.0;
		start = clock();
		for (int i = 0; i < reps; ++i) {
			C = A / B;
		}
		assert(C(0,0)==1);
		end = clock();
		assert( start != -1);
		assert( end != -1);
		printReportMsg("out-of-place matrix-matrix division", start, end, reps);


		// --------------------------------------------------- CHECK "CLOSE TO"
		A = 2.0;
		B = gMat2D<T>::rand(Nrows,Ncols)/static_cast<T>(1E11+2);
		start = clock();
		for (int i = 0; i < reps; ++i) {
			A.closeTo(B, (T)1E-10);
		}
		assert(A.closeTo(B, (T)1E-10));
		end = clock();
		printReportMsg("test of method closeto()", start, end, reps);

		// ------------------------------------------------------ TRANSPOSITION
		start = clock();
		for (int i = 0; i < reps; ++i) {
			gMat2D<T> D(Ncols, Nrows);
			A.transpose(D);
		}
		end = clock();
		printReportMsg("transposition", start, end, reps);
	}
	catch (gException& e) {

		cout << e.getMessage() << endl;
	}

}

void benchmarkGMATH(int Nrows, int Ncols, int reps) {
	clock_t start,end;

	gMat2D<T> *A, *B, *C;
	gVec<T> *X, *Y, *U;

	T value;

	try {
		X = new gVec<T>(Nrows);
		Y = new gVec<T>(Nrows);
		*X = 1;
		*X = 2;
		start = clock();
		for (int i = 0; i < reps; ++i) {
			value = dot(*X,*Y);
		}
		end = clock();
		printReportMsg("vector-vector product", start, end, reps);

		start = clock();
		A = new gMat2D<T>(Nrows, Ncols);

		U = new gVec<T>(Ncols);
		for (int i = 0; i < reps; ++i) {
			X = new gVec<T>(Nrows);
			dot(*A,*U,*X);
			delete X;
		}
		delete A;
		delete U;
		end = clock();
		printReportMsg("out-of-place matrix-vector product", start, end, reps);

		start = clock();
		A = new gMat2D<T>(Nrows, Ncols);
		B = new gMat2D<T>(Ncols, Nrows);
		for (int i = 0; i < reps; ++i) {
			C = new gMat2D<T>(Nrows, Nrows);
			dot(*A,*B,*C);
			delete C;
		}
		delete A;
		delete B;
		end = clock();
		printReportMsg("out-of-place matrix-matrix product", start, end, reps);
	}
	catch (gException& e) {

		cout << e.getMessage() << endl;
	}
}

void benchmarkNORM(int Nrows, int Ncols, int reps){
	clock_t start,end;

	//gMat2D<T> A(Nrows,Ncols), B(Ncols,Nrows), C(Nrows,Nrows);
	gVec<T> X(Nrows), U(Ncols);

	X = 1;
	U = 1;

	try {
		start = clock();
		for (int i = 0; i < reps; ++i) {
			norm(X, L0norm);
		}
		end = clock();
		printReportMsg("L0 norm", start, end, reps);

		start = clock();
		for (int i = 0; i < reps; ++i) {
			norm(X, L1norm);
		}
		end = clock();
		printReportMsg("L1 norm", start, end, reps);

		start = clock();
		for (int i = 0; i < reps; ++i) {
			norm(X, L2norm);
		}
		end = clock();
		printReportMsg("L2 norm", start, end, reps);
	}
	catch (gException& e) {

		cout << e.getMessage() << endl;
	}
}

void benchmarkKERNEL(int Nrows, int Ncols, int reps){
	// TO BE IMPLEMENTED
}

