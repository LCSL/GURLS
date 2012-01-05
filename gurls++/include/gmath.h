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


/**
\ingroup LinearAlgebra
\file
\brief Contains functions that implement some relevant BLAS level 1,2
or 3 and Lapack functionalities.

*/

#ifndef _GURLS_GMATH_H_
#define _GURLS_GMATH_H_

#include <cstring>
#include <cmath>

#include "gmat2d.h"
#include "gvec.h"
#include "exceptions.h"

extern "C" {
#include <cblas.h>
}

namespace gurls	{

template <typename T>
class gMat2D;

template <typename T>
class gVec;

//! Implements the standard GEMM routine from Level3 BLAS
/*! \fn void dot(const gMat2D<T>& A, const gMat2D<T>& B, gMat2D<T>& C);

General Matrix-Matrix multiplication of two single/double precision
real matrices A and B (the corresponding Matlab code is: C = A*B;).

 \param[in] A Left matrix
 \param[in] B Right matrix
 \param[out] C Product Matrix

*/
template <typename T>
void dot(const gMat2D<T>& A, const gMat2D<T>& B, gMat2D<T>& C);

//!  Implements the standard DOT routine from Level 1 BLAS
/*! \fn void dot(const gMat2D<T>& A, const gVec<T>& x, gVec<T>& y);

General Matrix-Vector multiplication of a single/double precision
matrix A with a vector x (the corresponding Matlab code is y = A*x;).

\param[in] A Input matrix
\param[in] x Input vector
\param[out] y Output vector

*/
template <typename T>
void dot(const gMat2D<T>& A, const gVec<T>& x, gVec<T>& y);

//! Implements the standard scalar product between vectors
/*! \fn T dot(const gVec<T>& x, const gVec<T>& y);

General routine from Level1 BLAS: n <-- x^T * y
 General Vector-Vector multiplication for single/double precision real data
*/
template <typename T>
T dot(const gVec<T>& x, const gVec<T>& y);


// =========================================================
// ================= TO BE DISCUSSED =======================
// =========================================================
/* Replicate a vector n times along the columns (or along the rows if transpose==true)
   If x is a vector of length N, then the output is an N-by-N matrix whose columns (or rows) are copies of x.
*/
template <typename T>
gMat2D<T> repmat(const gVec<T>& x, unsigned long n, bool transpose = false){

	T dataA[n*x.size()];
	const T* datax = x.getData();
	ulong rows, cols;

	if(transpose){
		rows = n;
		cols = x.getSize();
		for(int i = 0; i < n; i++) {
			memcpy(dataA+i*cols, datax, cols);
		}
	}else {
		rows = x.size();
		cols = n;
		for(int i = 0; i<rows; i++) {
			for(int j = 0; j < cols; j++) {
				*(dataA+i*cols+j) = *(datax+i);
			}
		}
	}

	gMat2D<T> A(dataA, rows, cols, true);
	return A;
}

enum InversionAlgorithm {LU,GaussJ};

/*
 Implements the LU decomposition usig LAPACK routines
 */
template <typename T>
void lu(gMat2D<T>& A);


/*
 Implements the LU factorization of a general M-by-N matrix A using partial
 pivoting with row interchanges,, using LAPACK routines.

 The factorization has the form  A = P * L * U  where  P  is a permutation
 matrix, L is lower triangular with unit diagonal elements (lower
 trapezoidal if m > n), and U is upper triangular (upper trapezoidal
 if m < n). This is the right-looking Level 3 BLAS version of the algorithm.
*/
template <typename T>
void lu(gMat2D<T>& A, gVec<int>& pv);

/*
 Implements the SVD decomposition of a general rectangular matrix: A = U*W*Vt
 A Matrix to be decomposed
 U Matrix of the left singular vectors
 W Vector containing the singular values of the decomposition
 Vt transposed matrix of the right singular vectors
*/
template <typename T>
void svd(const gMat2D<T>& A, gMat2D<T>& U, gVec<T>& W, gMat2D<T>& Vt);

/* Implements the computation of the eigenvalues of A using the LAPACK routine
 SGEEV with default computation of the right eigenvectors.
 */
template <typename T>
void eig(const gMat2D<T>& A, gMat2D<T>& V, gVec<T>& Wr, gVec<T>& Wi);

/* Implements the computation of the eigenvalues of A using the LAPACK routine
 SGEEV with default computation of the right eigenvectors.
 */
template <typename T>
void eig(const gMat2D<T>& A, gMat2D<T>& V, gVec<T>& W);


/*
 Implements the computation of the eigenvalues of A
*/
template <typename T>
void eig(const gMat2D<T>& A, gVec<T>& Wr, gVec<T>& Wi);

/*
Implements the computation of the eigenvalues of A
*/template <typename T>
void eig(const gMat2D<T>& A, gVec<T>& W);


/*
 Compute the inverse (A^-1) of a matrix A
*/
template <typename T>
void inv(const gMat2D<T>& A, gMat2D<T>& Ainv, InversionAlgorithm alg = LU);

/*
  Compute the pseudo-inverse of a non square matrix A
*/
template <typename T>
void pinv(const gMat2D<T>& A, gMat2D<T>& Ainv, T RCOND = 0);

/*
 Computes the Cholesky factorization of a symmetric,
 positive definite matrix using the LAPACK routine SPOTRF
*/
template <typename T>
void  cholesky(const gMat2D<T>& A, gMat2D<T>& L, bool upper = true);

}

#endif // _GURLS_GMATH_H_
