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

#include "gmath.h"

#include "gmat2d.h"
#include "gvec.h"

namespace gurls {

/**
  * Specialized version of dot for float vectors
  */
template <>
GURLS_EXPORT float dot(const gVec<float>& x, const gVec<float>& y)
{
    if ( x.getSize() != y.getSize() )
        throw gException(gurls::Exception_Inconsistent_Size);

    int n = x.getSize();
    int incr = 1;

    return sdot_(&n, const_cast<float*>(x.getData()), &incr, const_cast<float*>(y.getData()), &incr);
}

/**
  * Specialized version of dot for float vectors
  */
template <>
GURLS_EXPORT double dot(const gVec<double>& x, const gVec<double>& y)
{

    if ( x.getSize() != y.getSize() )
        throw gException(gurls::Exception_Inconsistent_Size);

    int n = x.getSize();
    int incr = 1;

    return ddot_(&n, const_cast<double*>(x.getData()), &incr, const_cast<double*>(y.getData()), &incr);
}


// ============ OUT-OF-PLACE MATRIX MULTIPLICATION ==================
/**
  * Specialized version of dot for float matrices
  */
template <>
GURLS_EXPORT void dot(const gMat2D<float>& A, const gMat2D<float>& B, gMat2D<float>& C)
{

    dot(A.getData(), B.getData(), C.getData(),
        A.rows(), A.cols(),
        B.rows(), B.cols(),
        C.rows(), C.cols(),
        CblasNoTrans, CblasNoTrans, CblasColMajor);
}

/**
  * Specialized version of dot for float matrices
  */
template <>
GURLS_EXPORT void dot(const gMat2D<double>& A, const gMat2D<double>& B, gMat2D<double>& C)
{

    dot(A.getData(), B.getData(), C.getData(),
        A.rows(), A.cols(),
        B.rows(), B.cols(),
        C.rows(), C.cols(),
        CblasNoTrans, CblasNoTrans, CblasColMajor);
}


// ============ OUT-OF-PLACE MATRIX-VECTOR MULTIPLICATION ==================

/**
  * Specialized version of dot for float matrices/vectors
  */
template <>
GURLS_EXPORT void dot(const gMat2D<float>& A, const gVec<float>& x, gVec<float>& y)
{
    if ( (A.cols() != x.getSize()) ||  (A.rows() != y.getSize()))
        throw gException(Exception_Inconsistent_Size);


    // y = alpha*A*x + beta*y
    float alpha = 1.0f;
    float beta = 0.0f;

    char transA = 'N';

    int m = A.rows();
    int n = A.cols();
    int lda = m;
    int inc = 1;

    sgemv_(&transA, &m, &n, &alpha, const_cast<float*>(A.getData()), &lda,
          const_cast<float*>(x.getData()), &inc, &beta, y.getData(), &inc);
}

/**
  * Specialized version of dot for float matrices/vectors
  */
template <>
GURLS_EXPORT void dot(const gMat2D<double>& A, const gVec<double>& x, gVec<double>& y){

    if ( (A.cols() != x.getSize()) ||  (A.rows() != y.getSize()))
        throw gException(Exception_Inconsistent_Size);


    // y = alpha*A*x + beta*y
    double alpha = 1.0;
    double beta = 0.0;

    char transA = 'N';

    int m = A.rows();
    int n = A.cols();
    int lda = m;
    int inc = 1;

    dgemv_(&transA, &m, &n, &alpha, const_cast<double*>(A.getData()), &lda,
          const_cast<double*>(x.getData()), &inc, &beta,
          y.getData(), &inc);

}

/**
  * Specialized version of lu for float matrices/vectors
  */
template <>
GURLS_EXPORT void lu(gMat2D<float>& A, gVec<int>& pv)
{
    unsigned int k = std::min(A.cols(), A.rows());

    if (pv.getSize() != k)
        throw gException("The lenghth of pv must be equal to the minimun dimension of A");

    int info;
    int m = A.rows();
    int n = A.cols();
    int lda = A.rows();

    sgetrf_(&m, &n, A.getData(), &lda, pv.getData(), &info);

    if(info <0)
        throw gException("LU factorization failed");
}

/**
  * Specialized version of lu for float matrices
  */
template <>
GURLS_EXPORT void lu(gMat2D<float>& A)
{
    gVec<int> pv(std::min(A.cols(), A.rows()));
    lu(A, pv);
}

/**
  * Specialized version of inv for float matrices
  */
template <>
GURLS_EXPORT void inv(const gMat2D<float>& A, gMat2D<float>& Ainv, InversionAlgorithm alg)
{
    Ainv = A;
    int k = std::min(Ainv.cols(), Ainv.rows());

    int info;
    int* ipiv = new int[k];

    int m = Ainv.rows();
    int n = Ainv.cols();
    int lda = Ainv.rows();

    sgetrf_(&m, &n, Ainv.getData(), &lda, ipiv, &info);

    float* work = new float[n];

    sgetri_(&m, Ainv.getData(), &lda, ipiv, work, &n, &info);

    delete[] ipiv;
    delete[] work;
}


template<>
GURLS_EXPORT float* pinv(const float* A, int rows, int cols, int& res_rows, int& res_cols, float* RCOND);

/**
  * Specialized version of pinv for float matrices
  */
template <>
GURLS_EXPORT void pinv(const gMat2D<float>& A, gMat2D<float>& Ainv, float RCOND)
{
    int r, c;
    float* inv = gurls::pinv<float>(A.getData(), A.rows(), A.cols(), r, c, &RCOND);

    Ainv.resize(r, c);
    gurls::copy(Ainv.getData(), inv, r*c);

    delete[] inv;
}

/**
  * Specialized version of svd for float matrices/vectors
  */
template <>
GURLS_EXPORT void svd(const gMat2D<float>& A, gMat2D<float>& U, gVec<float>& W, gMat2D<float>& Vt)
{
    float* Ubuf;
    float* Sbuf;
    float* Vtbuf;

    int Urows, Ucols;
    int Slen;
    int Vtrows, Vtcols;

    gurls::svd(A.getData(), Ubuf, Sbuf, Vtbuf,
             A.rows(), A.cols(),
             Urows, Ucols, Slen, Vtrows, Vtcols);


    U.resize(Urows, Ucols);
    copy(U.getData(), Ubuf, U.getSize());

    W.resize(Slen);
    copy(W.getData(), Sbuf, Slen);

    Vt.resize(Vtrows, Vtcols);
    copy(Vt.getData(), Vtbuf, Vt.getSize());

    delete [] Ubuf;
    delete [] Sbuf;
    delete [] Vtbuf;
}

/**
  * Specialized version of eig for float matrices/vectors
  */
template <>
GURLS_EXPORT void eig(const gMat2D<float>& A, gMat2D<float>& V, gVec<float>& Wr, gVec<float>& Wi)
{
    if (A.cols() != A.rows())
        throw gException("The input matrix A must be squared");

    float* Atmp = new float[A.getSize()];
    copy(Atmp, A.getData(), A.getSize());

    char jobvl = 'N', jobvr = 'V';
    int n = A.cols(), lda = A.cols(), ldvl = 1, ldvr = A.cols();
    int info, lwork = 4*n;
    float* work = new float[lwork];

    sgeev_(&jobvl, &jobvr, &n, Atmp, &lda, Wr.getData(), Wi.getData(), NULL, &ldvl, V.getData(), &ldvr, work, &lwork, &info);

    delete[] Atmp;
    delete[] work;

    if(info != 0)
    {
        std::stringstream str;
        str << "Eigenvalues/eigenVectors computation failed, error code " << info << ";" << std::endl;
        throw gException(str.str());
    }
}

/**
  * Specialized version of eig for float matrices/vectors
  */
template <>
GURLS_EXPORT void eig(const gMat2D<float>& A, gMat2D<float>& V, gVec<float>& W)
{
    gVec<float> tmp(W.getSize());
    tmp = 0;

    eig(A, V, W, tmp);
}

/**
  * Specialized version of eig for float matrices/vectors
  */
template <>
GURLS_EXPORT void eig(const gMat2D<float>& A, gVec<float>& Wr, gVec<float>& Wi)
{
    if (A.cols() != A.rows())
        throw gException("The input matrix A must be squared");

    float* Atmp = new float[A.getSize()];
    copy(Atmp, A.getData(), A.getSize());

    char jobvl = 'N', jobvr = 'N';
    int n = A.cols(), lda = A.cols(), ldvl = 1, ldvr = 1;
    int info, lwork = 4*n;
    float* work = new float[lwork];

    sgeev_(&jobvl, &jobvr, &n, Atmp, &lda, Wr.getData(), Wi.getData(), NULL, &ldvl, NULL, &ldvr, work, &lwork, &info);

    delete[] Atmp;
    delete[] work;

    if(info != 0)
    {
        std::stringstream str;
        str << "Eigenvalues/eigenVectors computation failed, error code " << info << ";" << std::endl;
        throw gException(str.str());
    }
}

/**
  * Specialized version of eig for float matrices/vectors
  */
template <>
GURLS_EXPORT void eig(const gMat2D<float>& A, gVec<float>& W)
{
    gVec<float> tmp = W;
    eig(A, W, tmp);
}


/**
  * Specialized version of cholesky for float matrices
  */
template <>
GURLS_EXPORT void cholesky(const gMat2D<float>& A, gMat2D<float>& L, bool upper)
{
    cholesky<float>(A.getData(), A.rows(), A.cols(), L.getData(), upper);
}

/**
  * Specialized version of set for float buffers
  */
template<>
void GURLS_EXPORT set(float* buffer, const float value, const int size, const int incr)
{
    int incx = 0;

    scopy_(const_cast<int*>(&size), const_cast<float*>(&value), &incx, buffer, const_cast<int*>(&incr));
}

/**
  * Specialized version of set for float buffers
  */
template<>
void GURLS_EXPORT set(float* buffer, const float value, const int size)
{
    set<float>(buffer, value, size, 1);
}

/**
  * Specialized version of set for double buffers
  */
template<>
void GURLS_EXPORT set(double* buffer, const double value, const int size, const int incr)
{
    int incx = 0;

    dcopy_(const_cast<int*>(&size), const_cast<double*>(&value), &incx, buffer, const_cast<int*>(&incr));

}

/**
  * Specialized version of set for double buffers
  */
template<>
void GURLS_EXPORT set(double* buffer, const double value, const int size)
{
    set<double>(buffer, value, size, 1);
}

/**
  * Specialized version of copy for float buffers
  */
template<>
void GURLS_EXPORT copy(float* dst, const float* src, const int size, const int dstIncr, const int srcIncr)
{
    scopy_(const_cast<int*>(&size), const_cast<float*>(src), const_cast<int*>(&srcIncr), dst, const_cast<int*>(&dstIncr));
}

/**
  * Specialized version of copy for float buffers
  */
template<>
void GURLS_EXPORT copy(float* dst, const float* src, const int size)
{
    int incr = 1;

    scopy_(const_cast<int*>(&size), const_cast<float*>(src), &incr, dst, &incr);
}

/**
  * Specialized version of copy for double buffers
  */
template<>
void GURLS_EXPORT copy(double* dst, const double* src, const int size, const int dstIncr, const int srcIncr)
{
    dcopy_(const_cast<int*>(&size), const_cast<double*>(src), const_cast<int*>(&srcIncr), dst, const_cast<int*>(&dstIncr));
}

/**
  * Specialized version of copy for double buffers
  */
template<>
void GURLS_EXPORT copy(double* dst, const double* src, const int size)
{
    int incr = 1;

    dcopy_(const_cast<int*>(&size), const_cast<double*>(src), &incr, dst, &incr);
}

/**
  * Specialized version of pinv for float buffers
  */
template<>
GURLS_EXPORT float* pinv(const float* A, int rows, int cols, int& res_rows, int& res_cols, float* RCOND)
{
    int M = rows;
    int N = cols;

    float* a = new float[rows*cols];
    copy<float>(a, A, rows*cols);

    int LDA = M;
    int LDB = std::max(M, N);
    int NRHS = LDB;


    // float* b = eye(LDB).getData()

    const int b_size = LDB*NRHS;
    float *b = new float[LDB*NRHS];
    set<float>(b, 0.f, b_size);
    set<float>(b, 1.f, std::min(LDB, NRHS), NRHS+1);

    float* S = new float[std::min(M,N)];
//    float condnum = 0.f; // The condition number of A in the 2-norm = S(1)/S(min(m,n)).

    float rcond = (RCOND == NULL)? (std::max(rows, cols)*FLT_EPSILON): *RCOND;

//    if (RCOND < 0)
//        RCOND = 0.f;

    int RANK = -1; // std::min(M,N);
    int LWORK = -1; //2 * (3*LDB + std::max( 2*std::min(M,N), LDB));
    float* WORK = new float[1];

    /*

    subroutine SGELSS 	( 	INTEGER  	M,
      INTEGER  	N,
      INTEGER  	NRHS,
      REAL,dimension( lda, * )  	A,
      INTEGER  	LDA,
      REAL,dimension( ldb, * )  	B,
      INTEGER  	LDB,
      REAL,dimension( * )  	S,
      REAL  	RCOND,
      INTEGER  	RANK,
      REAL,dimension( * )  	WORK,
      INTEGER  	LWORK,
      INTEGER  	INFO
     )

    */

    /*
   INFO:
   = 0:	successful exit
   < 0:	if INFO = -i, the i-th argument had an illegal value.
   > 0:	the algorithm for computing the SVD failed to converge;
     if INFO = i, i off-diagonal elements of an intermediate
     bidiagonal form did not converge to zero.
   */
    int INFO;

    /* Query and allocate the optimal workspace */
    sgelss_( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &rcond, &RANK, WORK, &LWORK, &INFO);
    LWORK = static_cast<int>(WORK[0]);
    delete [] WORK;
    WORK = new float[LWORK];

    sgelss_( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &rcond, &RANK, WORK, &LWORK, &INFO);

    // TODO: check INFO on exit
    //condnum = S[0]/(S[std::min(M, N)]-1);

    if(INFO != 0)
    {
        std::stringstream str;
        str << "Pinv failed, error code " << INFO << ";" << std::endl;
        throw gException(str.str());
    }

    delete [] S;
    delete [] WORK;
    delete [] a;

    res_rows = LDB;
    res_cols = NRHS;
    return b;

}

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
  * Specialized version of eq for double values
  */
template<>
GURLS_EXPORT bool eq(double val1, double val2)
{
    return (val1 >= val2-DBL_EPSILON && val1 <= val2+DBL_EPSILON );
}

/**
  * Specialized version of eq for float values
  */
template<>
GURLS_EXPORT bool eq(float val1, float val2)
{
    return ( val1 >= val2-FLT_EPSILON && val1 <= val2+FLT_EPSILON );
}

/**
  * Specialized version of gt for float values
  */
template<>
GURLS_EXPORT bool gt(double a, double b)
{
    return ((a - b) > ( std::min(fabs(a), fabs(b))* std::numeric_limits<double>::epsilon()));
}

/**
  * Specialized version of gt for double values
  */
template<>
GURLS_EXPORT bool gt(float a, float b)
{
    return ((a - b) > ( std::min(fabs(a), fabs(b))* std::numeric_limits<float>::epsilon()));
}

/**
  * Specialized version of lt for float values
  */
template<>
GURLS_EXPORT bool lt(double a, double b)
{
    return ((b - a) > ( std::max(fabs(a), fabs(b))* std::numeric_limits<double>::epsilon()));
}

/**
  * Specialized version of lt for double values
  */
template<>
GURLS_EXPORT bool lt(float a, float b)
{
    return ((b - a) > ( std::max(fabs(a), fabs(b))* std::numeric_limits<float>::epsilon()));
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

}
