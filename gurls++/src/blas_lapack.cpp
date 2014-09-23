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

#include "gurls++/blas_lapack.h"
#include "gurls++/exports.h"

namespace gurls {

/**
  * Specialized version of gemm for float buffers
  */
template<>
GURLS_EXPORT void gemm(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K, const float alpha, const float *A, const int lda,
          const float *B, const int ldb, const float beta, float *C, const int ldc)
{
    char transA = BlasUtils::charValue(TransA);
    char transB = BlasUtils::charValue(TransB);

    sgemm_(&transA, &transB,
          const_cast<int*>(&M), const_cast<int*>(&N), const_cast<int*>(&K),
          const_cast<float*>(&alpha), const_cast<float*>(A), const_cast<int*>(&lda),
          const_cast<float*>(B), const_cast<int*>(&ldb), const_cast<float*>(&beta),
          const_cast<float*>(C), const_cast<int*>(&ldc));

}

/**
  * Specialized version of gemm for double buffers
  */
template<>
GURLS_EXPORT void gemm(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K, const double alpha, const double *A, const int lda,
          const double *B, const int ldb, const double beta, double *C, const int ldc)
{
    char transA = BlasUtils::charValue(TransA);
    char transB = BlasUtils::charValue(TransB);

    dgemm_(&transA, &transB,
          const_cast<int*>(&M), const_cast<int*>(&N), const_cast<int*>(&K),
          const_cast<double*>(&alpha), const_cast<double*>(A), const_cast<int*>(&lda),
          const_cast<double*>(B), const_cast<int*>(&ldb), const_cast<double*>(&beta),
          const_cast<double*>(C), const_cast<int*>(&ldc));
}

/**
  * Specialized version of potrf_ for float buffers
  */
template<>
GURLS_EXPORT int potrf_(char *UPLO, int *n, float *a, int *lda , int *info)
{
    return spotrf_(UPLO, n, a, lda, info);
}

/**
  * Specialized version of potrf_ for double buffers
  */
template<>
GURLS_EXPORT int potrf_(char *UPLO, int *n, double *a, int *lda , int *info)
{
    return dpotrf_(UPLO, n, a, lda, info);
}

/**
  * Specialized version of axpy for float buffers
  */
template<>
GURLS_EXPORT void axpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
{
    saxpy_(const_cast<int*>(&N), const_cast<float*>(&alpha), const_cast<float*>(X), const_cast<int*>(&incX), Y, const_cast<int*>(&incY));
}

/**
  * Specialized version of axpy for double buffers
  */
template<>
GURLS_EXPORT void axpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY)
{
    daxpy_(const_cast<int*>(&N), const_cast<double*>(&alpha), const_cast<double*>(X), const_cast<int*>(&incX), Y, const_cast<int*>(&incY));
}

/**
  * Specialized version of dot for float buffers
  */
template <>
GURLS_EXPORT float dot(const int N, const float *X, const int incX, const float *Y, const int incY)
{
    return sdot_(const_cast<int*>(&N), const_cast<float*>(X), const_cast<int*>(&incX), const_cast<float*>(Y), const_cast<int*>(&incY));
}

/**
  * Specialized version of dot for double buffers
  */
template <>
GURLS_EXPORT double dot(const int N, const double *X, const int incX, const double *Y, const int incY)
{
    return ddot_(const_cast<int*>(&N), const_cast<double*>(X), const_cast<int*>(&incX), const_cast<double*>(Y), const_cast<int*>(&incY));
}

/**
  * Specialized version of nrm2 for float buffers
  */
template<>
GURLS_EXPORT float nrm2(const int N, const float* X, const int incX)
{
    return snrm2_(const_cast<int*>(&N), const_cast<float*>(X), const_cast<int*>(&incX));
}

/**
  * Specialized version of nrm2 for double buffers
  */
template<>
GURLS_EXPORT double nrm2(const int N, const double* X, const int incX)
{
    return dnrm2_(const_cast<int*>(&N), const_cast<double*>(X), const_cast<int*>(&incX));
}

/**
  * Specialized version of scal for float buffers
  */
template<>
GURLS_EXPORT void scal(const int N, const float alpha, float *X, const int incX)
{
    sscal_(const_cast<int*>(&N), const_cast<float*>(&alpha), X, const_cast<int*>(&incX));
}

/**
  * Specialized version of scal for double buffers
  */
template<>
GURLS_EXPORT void scal(const int N, const double alpha, double *X, const int incX)
{
    dscal_(const_cast<int*>(&N), const_cast<double*>(&alpha), X, const_cast<int*>(&incX));
}

/**
  * Specialized version of gemv for float buffers
  */
template<>
GURLS_EXPORT void gemv(const CBLAS_TRANSPOSE TransA,
          const int M, const int N, const float alpha, const float *A,
          const int lda, const float *X, const int incX,
          const float beta, float *Y, const int incY)
{
    char transA = BlasUtils::charValue(TransA);

    sgemv_(&transA, const_cast<int*>(&M), const_cast<int*>(&N),
          const_cast<float*>(&alpha), const_cast<float*>(A), const_cast<int*>(&lda),
          const_cast<float*>(X), const_cast<int*>(&incX), const_cast<float*>(&beta),
          const_cast<float*>(Y), const_cast<int*>(&incY));

}

/**
  * Specialized version of gemv for double buffers
  */
template<>
GURLS_EXPORT void gemv(const CBLAS_TRANSPOSE TransA,
          const int M, const int N, const double alpha, const double *A,
          const int lda, const double *X, const int incX,
          const double beta, double *Y, const int incY)
{
    char transA = BlasUtils::charValue(TransA);

    dgemv_(&transA, const_cast<int*>(&M), const_cast<int*>(&N),
          const_cast<double*>(&alpha), const_cast<double*>(A), const_cast<int*>(&lda),
          const_cast<double*>(X), const_cast<int*>(&incX), const_cast<double*>(&beta),
          const_cast<double*>(Y), const_cast<int*>(&incY));
}
/**
  * Specialized version of rot for float buffers
  */
template<>
GURLS_EXPORT void rot(int *N, float *X, int *incX, float *Y, int *incY, float *c, float *s)
{
srot_(N, X, incX, Y, incY, c, s);
}

/**
  * Specialized version of rot for double buffers
  */
template<>
GURLS_EXPORT void rot(int *N, double *X, int *incX, double *Y, int *incY, double *c, double *s)
{
drot_(N, X, incX, Y, incY, c, s);
}

/**
  * Specialized version of rotg for float buffers
  */
template<>
GURLS_EXPORT void rotg(float *a, float *b, float *c, float *s)
{
srotg_(a, b, c, s);
}

/**
  * Specialized version of rotg for double buffers
  */
template<>
GURLS_EXPORT void rotg(double *a, double *b, double *c, double *s)
{
drotg_(a, b, c, s);
}


/**
  * Specialized version of syev for float buffers
  */
template<>
GURLS_EXPORT void syev( char* jobz, char* uplo, int* n, float* a, int* lda, float* w, float* work, int* lwork, int* info)
{
    ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

/**
  * Specialized version of syev for double buffers
  */
template<>
GURLS_EXPORT void syev( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info)
{
    dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

/**
  * Specialized version of trsm for float buffers
  */
template <>
GURLS_EXPORT void trsm(const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N, const float alpha, const float *A, const int lda, float *B, const int ldb)
{
    char side = BlasUtils::charValue(Side);
    char uplo = BlasUtils::charValue(Uplo);
    char transA = BlasUtils::charValue(TransA);
    char diag = BlasUtils::charValue(Diag);

    strsm_(&side, &uplo, &transA, &diag, const_cast<int*>(&M), const_cast<int*>(&N), const_cast<float*>(&alpha), const_cast<float*>(A),
          const_cast<int*>(&lda), const_cast<float*>(B), const_cast<int*>(&ldb));

}

/**
  * Specialized version of trsm for double buffers
  */
template <>
GURLS_EXPORT void trsm(const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N, const double alpha, const double *A, const int lda, double *B, const int ldb)
{
    char side = BlasUtils::charValue(Side);
    char uplo = BlasUtils::charValue(Uplo);
    char transA = BlasUtils::charValue(TransA);
    char diag = BlasUtils::charValue(Diag);

    dtrsm_(&side, &uplo, &transA, &diag, const_cast<int*>(&M), const_cast<int*>(&N), const_cast<double*>(&alpha), const_cast<double*>(A),
          const_cast<int*>(&lda), const_cast<double*>(B), const_cast<int*>(&ldb));
}

/**
  * Specialized version of gesvd_ for float buffers
  */
template <>
GURLS_EXPORT int gesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info)
{
    return sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}

/**
  * Specialized version of gesvd_ for double buffers
  */
template <>
GURLS_EXPORT int gesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info)
{
    return dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}

/**
  * Specialized version of geqp3 for float buffers
  */
template<>
GURLS_EXPORT void geqp3( int *m, int *n, float *A, int *lda, int *jpvt, float *tau, float *work, int *lwork, int *info)
{
    sgeqp3_(m, n, A, lda, jpvt, tau, work, lwork, info);
}

/**
  * Specialized version of geqp3 for double buffers
  */
template<>
GURLS_EXPORT void geqp3( int *m, int *n, double *A, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info)
{
    dgeqp3_(m, n, A, lda, jpvt, tau, work, lwork, info);
}

/**
  * Specialized version of orogqr for float buffers
  */
template<>
GURLS_EXPORT void orgqr(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info)
{
    sorgqr_(m, n, k, a, lda, tau, work, lwork, info);
}

/**
  * Specialized version of orgqr for double buffers
  */
template<>
GURLS_EXPORT void orgqr(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info)
{
    dorgqr_(m, n, k, a, lda, tau, work, lwork, info);
}

/**
  * Specialized version of gelss for float buffers
  */
template<>
GURLS_EXPORT int gelss( int *m, int *n, int* nrhs, float *a, int *lda, float* b, int *ldb, float *s, float *rcond, int *rank, float *work, int *lwork, int *info)
{
    return sgelss_( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}

/**
  * Specialized version of gelss for double buffers
  */
template<>
GURLS_EXPORT int gelss( int *m, int *n, int* nrhs, double *a, int *lda, double* b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *info)
{
    return dgelss_( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}

/**
  * Specialized version of swap for float buffers
  */
template<>
GURLS_EXPORT void swap( int n, float *x, int incx, float *y, int incy)
{
    sswap_(&n, x, &incx, y, &incy);
}

/**
  * Specialized version of swap for double buffers
  */
template<>
GURLS_EXPORT void swap( int n, double *x, int incx, double *y, int incy)
{
    dswap_(&n, x, &incx, y, &incy);
}

}
