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


#ifndef BLAS_LAPACK_HPP
#define BLAS_LAPACK_HPP

/**
 * \ingroup LinearAlgebra
 * \file
 * \brief Contains template functions for BLAS level 1, 2, 3 and Lapack
 * routines.
 */

namespace gurls {

/**
  * Template function to call LAPACK *GELSS routines
  */
template<typename T>
int gelss( int *m, int *n, int* nrhs, T *a, int *lda, T* b, int *ldb, T *s, T *rcond, int *rank, T *work, int *lwork, int *info);

/**
  * Template function to call BLAS *SWAP routines
  */
template<typename T>
void swap( int n, T *x, int incx, T *y, int incy);

/**
  * Template function to call BLAS *GESVD routines
  */
template<typename T>
int gesvd_(char *jobu, char *jobvt, int *m, int *n, T *a, int *lda, T *s, T *u, int *ldu, T *vt, int *ldvt, T *work, int *lwork, int *info);

/**
  * Template function to call BLAS *SCAL routines
  */
template<typename T>
void scal(const int N, const T alpha, T *X, const int incX);

/**
  * Template function to call LAPACK *POTRF routines
  */
template<typename T>
int potrf_(char *UPLO, int *n, T *a, int *lda , int *info);

/**
  * Template function to call BLAS *AXPY routines
  */
template <typename T>
void axpy(const int N, const T alpha, const T *X, const int incX, T *Y, const int incY);

/**
  * Template function to call BLAS *DOT routines
  */
template <typename T>
T dot(const int N, const T *X, const int incX, const T *Y, const int incY);

/**
  * Template function to call BLAS *TRSM routines
  */
template <typename T>
void trsm(const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N, const T alpha, const T *A, const int lda, T *B, const int ldb);

/**
  * Template function to call BLAS *NRM2 routines
  */
template<typename T>
T nrm2(const int N, const T* X, const int incX);

/**
  * Template function to call BLAS *GEMV routines
  */
template<typename T>
void gemv(const CBLAS_TRANSPOSE TransA,
          const int M, const int N, const T alpha, const T *A, const int lda,
          const T *X, const int incX,
          const T beta, T *Y, const int incY);
		  
/**
  * Template function to call BLAS *ROT routines
  */
template<typename T>
void rot(int *N, T *X, int *incX, T *Y, int *incY, T *c, T *s);
		  
/**
  * Template function to call BLAS *ROTG routines
  */
template<typename T>
void rotg(T *a, T *b, T *c, T *s);

/**
  * Template function to call LAPACK *SYEV routines
  */
template<typename T>
void syev( char* jobz, char* uplo, int* n, T* a, int* lda, T* w, T* work, int* lwork, int* info);


/**
  * Template function to call BLAS *GEMM routines
  */
template<typename T>
void gemm(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K, const T alpha, const T *A, const int lda,
          const T *B, const int ldb,
          const T beta, T *C, const int ldc);

/**
  * Template function to call LAPACK *GEQP3 routines
  */
template<typename T>
void geqp3( int *m, int *n, T *A, int *lda, int *jpvt, T *tau, T *work, int *lwork, int *info);

/**
  * Template function to call LAPACK *ORGQR routines
  */
template<typename T>
void orgqr(int *m, int *n, int *k, T *a, int *lda, T *tau, T *work, int *lwork, int *info);

}

#endif
