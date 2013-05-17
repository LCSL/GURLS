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


/**
 * \ingroup LinearAlgebra
 * \file
 * \brief Contains functions that implement some relevant BLAS level 1,2
 * or 3 and Lapack functionalities.
 */

#ifndef _GURLS_GMATH_H_
#define _GURLS_GMATH_H_

#include <cstring>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <limits>

#include "exports.h"
#include "exceptions.h"

#ifdef _WIN32
#pragma warning(disable : 4290)
#endif

#include "blas_lapack.h"

namespace gurls	{

template <typename T>
class gMat2D;

template <typename T>
class gVec;

/**
  * \brief Implements the standard GEMM routine from Level3 BLAS
  * General Matrix-Matrix multiplication of two single/double precision
  * real matrices A and B (the corresponding Matlab code is: C = A*B;).
  *
  * \param[in] A Left matrix
  * \param[in] B Right matrix
  * \param[out] C Product Matrix
*/
template <typename T>
void dot(const gMat2D<T>& A, const gMat2D<T>& B, gMat2D<T>& C);

/**
  * Implements the standard DOT routine from Level 1 BLAS
  * General Matrix-Vector multiplication of a single/double precision
  * matrix A with a vector x (the corresponding Matlab code is y = A*x;).
  *
  * \param[in] A Input matrix
  * \param[in] x Input vector
  * \param[out] y Output vector
*/
template <typename T>
void dot(const gMat2D<T>& A, const gVec<T>& x, gVec<T>& y);

/**
  * Implements the standard scalar product between vectors
  * General routine from Level1 BLAS: n <-- x^T * y
  * General Vector-Vector multiplication for single/double precision real data
  */
template <typename T>
T dot(const gVec<T>& x, const gVec<T>& y);

/**
  * Replicates a vector n times along the columns (or along the rows if transpose==true)
  * If x is a vector of length N, then the output is an N-by-N matrix whose columns (or rows) are copies of x.
  */
template <typename T>
gMat2D<T> repmat(const gVec<T>& x, unsigned long n, bool transpose = false)
{
    gMat2D<T> A;

    const T* dataX = x.getData();
    const unsigned long len = x.getSize();

    if(transpose)
    {
        A.resize(n, len);

        for(int i = 0; i < n; ++i)
            copy(A.getData() + i, dataX, len, n, 1);
    }
    else
    {
        A.resize(len, n);

        for(int i = 0; i < n; ++i)
            copy(A.getData() + (i*len), dataX, len);
    }

    return A;
}

/**
  * \enum InversionAlgorithm
  * Implemented inversion algorithms
  */
enum InversionAlgorithm {LU,GaussJ};

/**
  * Implements the LU decomposition usig LAPACK routines
  */
template <typename T>
void lu(gMat2D<T>& A);


/**
  * Implements the LU factorization of a general M-by-N matrix A using partial
  * pivoting with row interchanges, using LAPACK routines.
  * The factorization has the form  A = P * L * U  where  P  is a permutation
  * matrix, L is lower triangular with unit diagonal elements (lower
  * trapezoidal if m > n), and U is upper triangular (upper trapezoidal
  * if m < n). This is the right-looking Level 3 BLAS version of the algorithm.
*/
template <typename T>
void lu(gMat2D<T>& A, gVec<int>& pv);

/**
  * Implements the SVD decomposition of a general rectangular matrix: A = U*W*Vt
  * \param A Matrix to be decomposed
  * \param U Matrix of the left singular vectors
  * \param W Vector containing the singular values of the decomposition
  * \param Vt transposed matrix of the right singular vectors
*/
template <typename T>
void svd(const gMat2D<T>& A, gMat2D<T>& U, gVec<T>& W, gMat2D<T>& Vt);

/**
  * Implements the computation of the eigenvalues of A using the LAPACK routine
  * SGEEV with default computation of the right eigenvectors.
 */
template <typename T>
void eig(const gMat2D<T>& A, gMat2D<T>& V, gVec<T>& Wr, gVec<T>& Wi);

/**
  * Implements the computation of the eigenvalues of A using the LAPACK routine
  * SGEEV with default computation of the right eigenvectors.
  */
template <typename T>
void eig(const gMat2D<T>& A, gMat2D<T>& V, gVec<T>& W);


/**
  * Implements the computation of the eigenvalues of A
  */
template <typename T>
void eig(const gMat2D<T>& A, gVec<T>& Wr, gVec<T>& Wi);

/**
  * Implements the computation of the eigenvalues of A
  */
template <typename T>
void eig(const gMat2D<T>& A, gVec<T>& W);


/**
  * Computes the inverse \f$A^-{1}\f$ of a matrix \f$A\f$
  */
template <typename T>
void inv(const gMat2D<T>& A, gMat2D<T>& Ainv, InversionAlgorithm alg = LU);

/**
  * Computes the pseudo-inverse of a non square matrix A
  */
template <typename T>
void pinv(const gMat2D<T>& A, gMat2D<T>& Ainv, T RCOND = 0);

/**
  * Computes the Cholesky factorization of a symmetric,
  * positive definite matrix using the LAPACK routine SPOTRF
  */
template <typename T>
void  cholesky(const gMat2D<T>& A, gMat2D<T>& L, bool upper = true);

/**
  * Sets elements of a vector to a specified value
  *
  * \param buffer vector to be set
  * \param value value to set
  * \param size number of copy operations to be performed
  * \param incr number of element to skip after every copy (if incr == 1 all elements will be set)
  */
template<typename T>
void set(T* buffer, const T value, const int size, const int incr);

/**
  * Sets all elements of a vector to a specified value
  *
  * \param buffer vector to be set
  * \param value value to set
  * \param size buffer length
  */
template<typename T>
void set(T* buffer, const T value, const int size)
{
    for(T *it = buffer, *end = buffer+size; it != end; ++it)
        *it = value;
}

/**
  * Specialized version of set for float buffers
  */
template<>
void GURLS_EXPORT set(float* buffer, const float value, const int size);

/**
  * Specialized version of set for double buffers
  */
template<>
void GURLS_EXPORT set(double* buffer, const double value, const int size);

/**
  * Copies element form one vector to another one
  *
  * \param dst destination vector
  * \param src source vector
  * \param size vectors length
  */
template<typename T>
void copy(T* dst, const T* src, const int size)
{
    memcpy(dst, src, size*sizeof(T));
}

/**
  * Specialized version of copy for float buffers
  */
template<>
GURLS_EXPORT void copy(float* dst, const float* src, const int size);

/**
  * Specialized version of copy for double buffers
  */
template<>
GURLS_EXPORT void copy(double* dst, const double* src, const int size);

/**
  * Copies element form one vector to another one
  *
  * \param dst destination vector
  * \param src source vector
  * \param size number of copies to be performed
  * \param dstIncr iteration increment on destination buffer
  * \param srcIncr iteration increment on source buffer
  */
template<typename T>
void copy(T* dst, const T* src, const int size, const int dstIncr, const int srcIncr)
{
    if(dstIncr == 1 && srcIncr == 1)
        copy(dst, src, size);

    else
    {
        const T* s_it = src;

        for(T *d_it = dst, *end = dst+(size*dstIncr); d_it != end; d_it += dstIncr, src+=srcIncr)
            *d_it = *s_it;
    }
}

/**
  * Specialized version of copy for float buffers
  */
template<>
GURLS_EXPORT void copy(float* dst, const float* src, const int size, const int dstIncr, const int srcIncr);

/**
  * Specialized version of copy for double buffers
  */
template<>
GURLS_EXPORT void copy(double* dst, const double* src, const int size, const int dstIncr, const int srcIncr);

/**
  * Generates a submatrix from an input matrix
  *
  * \param dst output submatrix
  * \param src input matrix
  * \param src_Rows number of rows of the input matrix
  * \param sizeRows number of rows of the output submatrix
  * \param sizeCols number of columns of the output submatrix
  * \param indices_rows vector containing the row indices to copy (length must be == sizeRows)
  * \param indices_cols vector containing the column indices to copy (length must be == sizeCols)
  */
template<typename T>
void copy_submatrix(T* dst, const T* src, const int src_Rows, const int sizeRows, const int sizeCols, unsigned long *indices_rows, unsigned long *indices_cols)
{
  int t=0;
  for (int j = 0; j < sizeCols; ++j)
    for (int i = 0; i < sizeRows; ++i)
    {
        dst[t] = src[ indices_rows[i] + indices_cols[j]*src_Rows ];
        ++t;
    }
}

/**
  * Generates a submatrix from an input matrix
  *
  * \param matrix input matrix
  * \param mRows number of rows of the input matrix
  * \param mCols number of columns of the input matrix
  * \param rowsIndices vector containing the row indices to copy
  * \param nIndices length of the indices vector
  * \param submat output submatrix
  */
template<typename T>
void subMatrixFromRows(const T* matrix, const int mRows, const int mCols, const unsigned long* rowsIndices, const int nIndices, T* submat)
{
    if(mRows < nIndices)
        throw gException(Exception_Inconsistent_Size);

    for(const unsigned long *it = rowsIndices, *end = rowsIndices+nIndices; it != end; ++it)
        copy(submat + (it-rowsIndices), matrix+(*it), mCols, nIndices, mRows);
}


/**
  * Generates a submatrix from an input matrix
  *
  * \param matrix input matrix
  * \param mRows number of rows of the input matrix
  * \param mCols number of columns of the input matrix
  * \param colsIndices vector containing the columns indices to copy
  * \param nIndices length of the indices vector
  * \param submat output submatrix
  */
template<typename T>
void subMatrixFromColumns(const T* matrix, const int mRows, const int mCols, const unsigned long* colsIndices, const int nIndices, T* submat)
{
    if(mCols < nIndices)
        throw gException(Exception_Inconsistent_Size);

    for(const unsigned long *it = colsIndices, *end = colsIndices+nIndices; it != end; ++it)
        copy(submat+(mRows*(it-colsIndices)), matrix+(mRows*(*it)), mRows);
}

/**
  * Zeroes on the lower triangle of a matrix
  *
  * \param matrix input matrix
  * \param rows number of rows
  * \param cols number of columns
  */
template <typename T>
void clearLowerTriangular(T* matrix, int rows, int cols)
{
    T* it = matrix;
    const T zero = static_cast<T>(0.0);

    for (int i = 1; i <= cols; ++i)
    {
        const int len = rows-i;

        if(len > 0)
        {
            it+=i;

            set(it, zero, len);

            it+=len;
        }
    }
}

/**
  * Computes the pseudo-inverse of a matrix
  *
  * \param A input matrix
  * \param rows number of rows
  * \param cols number of columns
  * \param res_rows on exit contains the number of rows of the pseudoinverse matrix
  * \param res_cols on exit contains the number of columns of the pseudoinverse matrix
  * \param RCOND used to determine the effective rank of A. If *RCOND < 0, machine precision is used
  * \return the pseudoinverse matrix
  */
template <typename T>
T* pinv(const T* A, const int rows, const int cols, int& res_rows, int& res_cols, T* RCOND = NULL)
{
    int M = rows;
    int N = cols;

    T* a = new T[rows*cols];
    copy(a, A, rows*cols);

    int LDA = M;
    int LDB = std::max(M, N);
    int NRHS = LDB;

    const int b_size = LDB*NRHS;
    T *b = new T[LDB*NRHS];

    set(b, (T)0.0, b_size);
    set(b, (T)1.0, std::min(LDB, NRHS), NRHS+1);

    T* S = new T[std::min(M,N)];
    T rcond = (RCOND == NULL)? (std::max(rows, cols)*std::numeric_limits<T>::epsilon()): *RCOND;

    int RANK = -1; // std::min(M,N);
    int LWORK = -1; //2 * (3*LDB + std::max( 2*std::min(M,N), LDB));
    T* WORK = new T[1];

    int INFO;

    /* Query and allocate the optimal workspace */
    gelss( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &rcond, &RANK, WORK, &LWORK, &INFO);
    LWORK = static_cast<int>(WORK[0]);
    delete [] WORK;
    WORK = new T[LWORK];

    gelss( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &rcond, &RANK, WORK, &LWORK, &INFO);

    delete [] S;
    delete [] WORK;
    delete [] a;

    if(INFO != 0)
    {
        delete [] b;

        std::stringstream str;
        str << "Pinv failed, error code " << INFO << ";" << std::endl;
        throw gException(str.str());
    }

    res_rows = LDB;
    res_cols = NRHS;

    return b;
}

/**
  * Transpose a matrix.
  *
  * \param matrix input matrix to be transposed
  * \param rows number of rows of the input matrix
  * \param cols number of columns of the input matrix
  * \param transposed on exit contains the transposed matrix. Must be already initialized
  * with a number of rows and columns equal to the number of columns and rows of the input matrix
  */
template <typename T>
void transpose(const T* matrix, const int rows, const int cols, T* transposed)
{
    T* d1 = transposed;
    const T* d0 = matrix;
    const int N = rows;

    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            *d1++=*(d0+c*N+r);
}

/**
  * General Matrix-Matrix multiplication of two single/double precision
  * real matrices A and B.
  *
  * \param[in] A Left matrix
  * \param[in] B Right matrix
  * \param[out] C Product Matrix
  * \param[in] A_rows numer of rows of the matrix A
  * \param[in] A_cols numer of columns of the matrix A
  * \param[in] B_rows numer of rows of the matrix B
  * \param[in] B_cols numer of columns of the matrix B
  * \param[in] C_rows numer of rows of the matrix C
  * \param[in] C_cols numer of columns of the matrix C
  * \param[in] TransA Transposition option form matrix A
  * \param[in] TransB Transposition option form matrix B
  * \param[in] Order order option for matrix multiplication
*/
template <typename T>
void dot(const T* A, const T* B, T* C,
         int A_rows, int A_cols,
         int B_rows, int B_cols,
         int C_rows, int C_cols,
         const CBLAS_TRANSPOSE TransA,
         const CBLAS_TRANSPOSE TransB,
         const CBLAS_ORDER Order)
{

    bool transposeA = (TransA == CblasTrans || TransA == CblasConjTrans);
    bool transposeB = (TransB == CblasTrans || TransB == CblasConjTrans);

    if(transposeA)
        std::swap(A_rows, A_cols);

    if(transposeB)
        std::swap(B_rows, B_cols);

    if ((C_rows != A_rows) || (C_cols != B_cols))
        throw gException(gurls::Exception_Inconsistent_Size);

    if (A_cols != B_rows)
        throw gException(gurls::Exception_Inconsistent_Size);

    const T alpha = 1.0;
    const T beta = 0.0;

    int lda, ldb, ldc;

    int k = A_cols;

    // mxk kxn mxn

    switch(Order)
    {
    case CblasColMajor:

        lda = transposeA? k: A_rows;
        ldb = transposeB? B_cols: k;
        ldc = C_rows;

        gemm(TransA, TransB, A_rows, B_cols, k, alpha, A, lda, B, ldb, beta, C, ldc);
        break;

    case CblasRowMajor:

        lda = transposeA? A_rows: k;
        ldb = transposeB? k: B_cols;
        ldc = C_cols;

        gemm(TransB, TransA, B_cols, A_rows, k, alpha, B, ldb, A, lda, beta, C, ldc);
        break;

    default:
        lda = ldb = ldc = 0;
        assert(0);
    }
}

/**
  * Computes the Cholesky factorization of a symmetric,
  * positive definite matrix using the LAPACK routine SPOTRF
  *
  * \param matrix input matrix
  * \param rows number of rows of the input matrix
  * \param cols number of columns of the input matrix
  * \param upper whether to store the upper or the lower triangle of the result matrix
  * \param resyult Cholesky factorization of the input matrix
  */
template<typename T>
void cholesky(const T* matrix, const int rows, const int cols, T* result, bool upper = true)
{
    copy(result, matrix, rows*cols);

    int LDA = rows;
    int nc = cols;
    char UPLO = upper? 'U':'L';
    int info;

    potrf_(&UPLO, &nc, result, &LDA, &info);

    if(info != 0)
    {
        std::stringstream str;
        str << "Cholesky factorization failed, error code " << info << ";" << std::endl;
        throw gException(str.str());
    }

    clearLowerTriangular(result, rows, cols);
}

/**
  * Template function to call LAPACK *POTRF routines
  */
template<typename T>
int potrf_(char *UPLO, int *n, T *a, int *lda , int *info);

/**
  * Computes the element-by-element multiplicative inverse of an input matrix
  *
  * \param matrix input matrix
  * \param len number of matrix elements
  */
template <typename T>
void setReciprocal(T* matrix, const int len)
{
    const T one = static_cast<T>(1.0);

    for(T *it = matrix, *end = matrix+len; it != end; ++it)
        *it = one / *it;
}

/**
  * Returns a squared matrix initialized in the diagonal with values from a vector
  *
  * \param vector input vector
  * \param len vector length
  * \param result On exit it contains a len-by-len matrix with vector on the diagonal
  */
template <typename T>
void diag(T* vector, const int len, T* result)
{
    const T zero = static_cast<T>(0.0);

    set(result, zero, len*len);
    copy(result, vector, len, len+1, 1);
}

/**
  * Multiplication of two scalars
  */
template <typename T>
T mul(T a, T b)
{
    return a*b;
}

/**
  * Division of two scalars
  */
template <typename T>
T div(T a, T b)
{
    return a/b;
}

/**
  * "Equals" operator between two scalars
  */
template<typename T>
bool eq(T val1, T val2)
{
    return val1 == val2;
}

/**
  * "Equals" operator between two scalars, specialized for float values
  */
template<>
GURLS_EXPORT bool eq(float val1, float val2);

/**
  * "Equals" operator between two scalars, specialized for double values
  */
template<>
GURLS_EXPORT bool eq(double val1, double val2);

/**
  * "Greater than" operator between two scalars
  */
template<typename T>
bool gt(T a, T b)
{
    return a > b;
}

/**
  * "Greater than" operator between two scalars, specialized for float values
  */
template<>
GURLS_EXPORT bool gt(float a, float b);

/**
  * "Greater than" operator between two scalars, specialized for double values
  */
template<>
GURLS_EXPORT bool gt(double a, double b);

/**
  * "Greater than or equals" operator between two scalars
  */
template<typename T>
bool gte(T a, T b)
{
    return gt(a,b) || eq(a,b);
}

/**
  * "Less or equal than" operator between two scalars
  */
template<typename T>
bool le(T a, T b)
{
    return !gt(a,b);
}

/**
  * "Less than" operator between two scalars
  */
template<typename T>
bool lt(T a, T b)
{
    return a < b;
}

/**
  * "Less than" operator between two scalars, specialized for float values
  */
template<>
GURLS_EXPORT bool lt(float a, float b);

/**
  * "Less than" operator between two scalars, specialized for double values
  */
template<>
GURLS_EXPORT bool lt(double a, double b);

/**
  * Applies an element by element binary operation \f$op\f$ over two vectors (\f$result = A op B\f$)
  *
  * \param A first input vector
  * \param B second input vector
  * \param result vector
  * \param len vectors length
  * \param op binary operation to perform
  */
template <typename T>
void binOperation(const T* A, const T* B, T* result, const int len, T(*op)(T,T))
{
    const T *a_it = A, *a_end = A+len;
    const T *b_it = B;
    T *r_it = result;

    while(a_it != a_end)
    {
        *r_it = (*op)((*a_it),(*b_it));

        ++a_it;
        ++b_it;
        ++r_it;
    }
}

/**
  * Element by element multiplication of two vectors
  *
  * \param A first input vector
  * \param B second input vector
  * \param result vector
  * \param len vectors length
  */
template <typename T>
void mult(const T* A, const T* B, T* result, const int len)
{
    binOperation<T>(A, B, result, len, &mul);
}

/**
  * Element by element division of two vectors
  *
  * \param A first input vector
  * \param B second input vector
  * \param result vector
  * \param len vectors length
  */
template <typename T>
void rdivide(const T* A, const T* B, T* result, const int len)
{
    binOperation<T>(A, B, result, len, &div);
}

/**
  * Sums all elements along the rows of a matrix
  *
  * \param A input matrix
  * \param result vecotr of length A_cols containing sums for each row of the matrix
  * \param A_rows number of rows of matrix A
  * \param A_cols number of columns of matrix A
  * \param res_length results vector length (MUST be == A_cols)
  */
template <typename T>
void sum(const T* A, T* result, const int A_rows, const int A_cols, const int res_length) throw (gException)
{
    if(A_cols != res_length)
        throw gException(Exception_Inconsistent_Size);

    const T *a_it = A, *a_end;

    const T zero = static_cast<T>(0.0);

    for(T *r_it = result, *r_end = result+A_cols; r_it != r_end; ++r_it)
    {
        *r_it = zero;
        a_end = a_it+A_rows;

        while(a_it != a_end)
        {
            *r_it += *a_it;

            ++a_it;
        }
    }
}

/**
  * Sums all elements along the columns of a matrix
  *
  * \param A input matrix
  * \param result vecotr of length A_rows containing sums for each column of the matrix
  * \param A_rows number of rows of matrix A
  * \param A_cols number of columns of matrix A
  */
template <typename T>
void sum_col(const T* A, T* result, const int A_rows, const int A_cols) throw (gException)
{
    const T *a_it = A, *a_end;

    const T zero = static_cast<T>(0.0);

    for(T *r_it = result, *r_end = result+A_rows; r_it != r_end; ++r_it)
    {
        *r_it = zero;

        for(a_it = A+(r_it-result), a_end = a_it+((A_cols-1)*A_rows); a_it != a_end; a_it += A_rows)
            *r_it += *a_it;
    }
}

/**
  * Computes the mean values along the rows of a matrix
  *
  * \param A input matrix
  * \param result vecotr of length A_cols containing mean values for each row of the matrix
  * \param A_rows number of rows of matrix A
  * \param A_cols number of columns of matrix A
  * \param res_length results vector length (MUST be == A_cols)
  */
template <typename T>
void mean(const T* A, T* result, const int A_rows, const int A_cols, const int res_length) throw (gException)
{
    sum(A, result, A_rows, A_cols, res_length);
    scal(res_length, (T)(1.0/A_rows), result, 1);
}

/**
  * Coputes the smallest elements along the rows of a matrix
  *
  * \param A input matrix
  * \param result vecotr of length A_cols containing the smallest elements for each row of the matrix
  * \param A_rows number of rows of matrix A
  * \param A_cols number of columns of matrix A
  * \param res_length results vector length (MUST be == A_cols)
  */
template <typename T>
void argmin(const T* A, unsigned long* result, const int A_rows, const int A_cols, const int res_length) throw (gException)
{
    if(A_cols != res_length)
        throw gException(Exception_Inconsistent_Size);

    const T *a_it = A;

    for(unsigned long *r_it = result, *r_end = result+A_cols; r_it != r_end; ++r_it, a_it += A_rows)
        *r_it = (std::min_element(a_it, a_it+A_rows) - a_it);
}

/**
 * Returns a subvector of an input vector, containing elements of the input vector whose index is contained into an indices vector
 *
 * \param locs indices vector
 * \param src input vector
 * \param locs_len length of the indices vector
 * \param src_len length of the input vector
 * \param result on exit it contains a subvector with all elements of the input vector whose index is contained into \c locs
 */
template <typename T>
void copyLocations(const unsigned long* locs, const T* src, const int locs_len, const int src_len, T* result)
{
    T* ptr_v = result;

    int val;
    for(const unsigned long* l_it = locs, *l_end=locs+locs_len; l_it != l_end; ++l_it, ++ptr_v)
    {
        val = *l_it;
        if((val < 0) || (val > src_len))
            throw gException(gurls::Exception_Index_Out_of_Bound);

        *ptr_v = src[val];
    }
}

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
  * Computes the sum of all elements of a vector
  *
  * \param V input vector
  * \param len vector length
  * \return sum of all elements of the input vector
  */
template <typename T>
T sumv(const T* V, const int len) throw (gException)
{
    T y = 1;
    return dot(len, V, 1, &y, 0);
}

/**
  * Computes RLS estimator given the singular value decomposition of the kernel matrix
  *
  * \param Q eigenvectors of the kernel matrix
  * \param L eigenvalues of the kernel matrix
  * \param Qty result of the matrix multiplication of the transpose of Q times the labels vector Y \f$(Q^T Y)\f$
  * \param C on exit contains the rls coefficients matrix
  * \param lambda regularization parameter
  * \param n number of training samples
  * \param Q_rows number of rows of the matrix Q
  * \param Q_cols number of columns of the matrix Q
  * \param L_length number of elements of the vector L
  * \param Qty_rows number of rows of the matrix Qty
  * \param Qty_cols number of columns of the matrix Qty
  */
template<typename T>
void rls_eigen(const T* Q, const T* L, const T* Qty, T* C, const T lambda, const int n,
             const int Q_rows, const int Q_cols,
             const int L_length,
             const int Qty_rows, const int Qty_cols)//  throw (gException)
{
    T* work = new T [(Q_rows+1)*L_length];

    rls_eigen(Q, L, Qty, C, lambda, n, Q_rows, Q_cols, L_length, Qty_rows, Qty_cols, work);

    delete [] work;
}

/**
  * Computes RLS estimator given the singular value decomposition of the kernel matrix
  *
  * \param Q eigenvectors of the kernel matrix
  * \param L eigenvalues of the kernel matrix
  * \param Qty result of the matrix multiplication of the transpose of Q times the labels vector Y \f$(Q^T Y)\f$
  * \param C on exit contains the rls coefficients matrix
  * \param lambda regularization parameter
  * \param n number of training samples
  * \param Q_rows number of rows of the matrix Q
  * \param Q_cols number of columns of the matrix Q
  * \param L_length number of elements of the vector L
  * \param Qty_rows number of rows of the matrix Qty
  * \param Qty_cols number of columns of the matrix Qty
  * \param work Work buffer of length L_length*(Q_rows+1)
  */
template<typename T>
void rls_eigen(const T* Q, const T* L, const T* Qty, T* C, const T lambda, const int n,
             const int Q_rows, const int Q_cols,
             const int L_length,
             const int Qty_rows, const int Qty_cols, T* work)//  throw (gException)
{
    //function C = rls_eigen(Q,L,QtY,lambda,n)

    //sQ = size(Q,1); -> Q_rows

    //L = L + n*lambda;
    T* L1 = work; // size L_length

    set(L1, n*lambda , L_length);

    axpy(L_length, (T)1.0, L, 1, L1, 1);

    //L = L.^(-1);
    setReciprocal(L1, L_length);

    // L = diag(L)
    T* QL = work+L_length; // size Q_rows*L_length

    copy(QL, Q, Q_rows*Q_cols);
    for(int i=0; i< Q_cols; ++i)
        scal(Q_rows, L1[i], QL+(Q_rows*i), 1);

    //C = (Q*L)*QtY;
    dot(QL, Qty, C, Q_rows, L_length, Qty_rows, Qty_cols, Q_rows, Qty_cols, CblasNoTrans, CblasNoTrans, CblasColMajor);
}

/**
  * Computes a "signum vector" of the same size as an input vector, where each element is:
  *  - 1 if the corresponding element of the input vector is greater than zero
  *  - 0 if the corresponding element of the input vector is zero
  *  - -1 if the corresponding element of the input vector is less than zero
  *
  * \param vector input vector
  * \param size number of elements of the input vector
  * \return the signum vector
  */
template<typename T>
T* sign(const T* vector, const int size)
{
    const T* v_it = vector;
    const T* v_end = vector+size;

    T* ret = new T[size];
    T* r_it = ret;

    const T zero = (T)0.0;

    while(v_it != v_end)
    {
        *r_it = static_cast<T>((eq(*v_it, zero)? 0.0 : ((*v_it > 0.0)? 1.0 : -1.0)));

        ++v_it;
        ++r_it;
    }

    return ret;
}

/**
  * Compares element by element two vectors using a binary predicate,
  * and returns a vector where each element is:
  * - 1 if comparison predicate was verified on the pair of elements at the corresponging index
  * - 0 if comparison predicate was not verified on the pair of elements at the corresponging index
  *
  * \param vector1 first input vector
  * \param vector2 second input vector
  * \param size number of elements of the input vectors
  * \param pred binary predicate used for comparison
  * \return the results vector
  */
template<typename T>
T* compare(const T* vector1, const T* vector2, const int size, bool(*pred)(T,T))
{
    const T* v1_end = vector1+size;

    T* ret = new T[size];
    T* r_it = ret;

    const T zero = static_cast<T>(0.0);
    const T one = static_cast<T>(1.0);

    for(const T *v1_it = vector1, *v2_it = vector2; v1_it != v1_end; ++v1_it, ++v2_it, ++r_it)
        *r_it = (*pred)(*v1_it, *v2_it)? one: zero;

    return ret;
}


/**
  * Compares each element of a vector with a threshold using a binary predicate
  * and returns a vector where each element is:
  * - 1 if comparison predicate was verified on the pair (element at the corresponging index, threshold)
  * - 0 if comparison predicate was not verified on the pair (element at the corresponging index, threshold)
  *
  * \param vector input vector
  * \param thr threshold
  * \param size number of elements of the input vector
  * \param pred binary predicate used for comparison
  * \return the results vector
  */
template<typename T>
T* compare(const T* vector, const T thr, const int size, bool(*pred)(T,T))
{
    const T* v1_end = vector+size;

    T* ret = new T[size];
    T* r_it = ret;

    const T zero = static_cast<T>(0.0);
    const T one = static_cast<T>(1.0);

    for(const T *v1_it = vector; v1_it != v1_end; ++v1_it, ++r_it)
        *r_it = (*pred)(*v1_it, thr)? one: zero;

    return ret;
}


/**
  * Sorts the elements of a matrix along the columns
  *
  * \param M Input matrix
  * \param rows Matrix rows
  * \param cols Matrix columns
  * \param pred Binary predicate used for comparison
  * \param values On ouptut contains the ordered values matrix
  * \param indices On ouptut contains the ordered indices matrix
  */
template<typename T>
void sort(const T* M, const unsigned long rows, const unsigned long cols, bool(*pred)(T,T), T* values, unsigned long* indices)
{
    typedef std::multimap<T, unsigned long, bool(*)(T,T) > MapType;

    for(unsigned long i=0; i< rows; ++i)
    {
        MapType data(pred);

        for (unsigned long j = 0; j < cols; ++j)
        {
            const unsigned long index = i+(rows*j);
            data.insert( std::pair<T,unsigned long>(M[index], j));
        }

        typename MapType::iterator it = data.begin();

        for (unsigned long j = 0; j < cols; ++j, ++it)
        {
            const unsigned long index = i+(rows*j);

            if(values != NULL)
                values[index] = it->first;

            if(indices != NULL)
                indices[index] = it->second;
        }
    }
}

/**
  * Returns the indices of the largest elements along different dimensions of a matrix.
  *
  * \param A input matrix
  * \param A_rows number of rows of the input matrix
  * \param A_cols number of columns of the input matrix
  * \param ind vector containing computed indices, length must be A_cols if dimension == 1,  or A_rows if dimension == 2
  * \param work work buffer of size >= 0 if dimension == 1, or size >= (A_rows*A_cols) if dimension == 2
  * \param dimension the dimension along which largest elements have to be computed
  * \return the results vector
  */
template<typename T>
void indicesOfMax(const T* A, const int A_rows, const int A_cols, unsigned long* ind, T* work, const int dimension) throw (gException)
{
    const T *m_it;
    int m_rows;
    int m_cols;

    switch(dimension)
    {
    case 1:
        m_it = A;
        m_rows = A_rows;
        m_cols = A_cols;
        break;
    case 2:
        transpose(A, A_rows, A_cols, work);
        m_it = work;
        m_rows = A_cols;
        m_cols = A_rows;
        break;
    default:
        throw gException(gurls::Exception_Illegal_Argument_Value);
    }

    for(unsigned long *r_it = ind, *r_end = ind+m_cols; r_it != r_end; ++r_it, m_it += m_rows)
        *r_it = (std::max_element(m_it, m_it+m_rows) - m_it);

}

/**
  * Returns the largest elements along different dimensions of a matrix.
  *
  * \param A input matrix
  * \param A_rows number of rows of the input matrix
  * \param A_cols number of columns of the input matrix
  * \param maxv vector containing computed values, length must be A_cols if dimension == 1,  or A_rows if dimension == 2
  * \param work work buffer of size >= 0 if dimension == 1, or size >= (A_rows*A_cols) if dimension == 2
  * \param dimension the dimension along which largest elements have to be computed
  * \return the results vector
  */
template<typename T>
void maxValues(const T* A, const int A_rows, const int A_cols, T* maxv, T* work, const int dimension) throw (gException)
{
    const T *m_it;
    int m_rows;
    int m_cols;

    switch(dimension)
    {
    case 1:
        m_it = A;
        m_rows = A_rows;
        m_cols = A_cols;
        break;
    case 2:
        transpose(A, A_rows, A_cols, work);
        m_it = work;
        m_rows = A_cols;
        m_cols = A_rows;
        break;
    default:
        throw gException(gurls::Exception_Illegal_Argument_Value);
    }

    for(T *r_it = maxv, *r_end = maxv+m_cols; r_it != r_end; ++r_it, m_it += m_rows)
        *r_it = *std::max_element(m_it, m_it+m_rows);

}

/**
  * Builds array of possible values for the regularization parameter,
  * generating a geometric series from the values in EIGVALS
  * Internal function, not to be called from gurls
  *
  * \param eigvals vector containing the eigenvalues of \f$X^T X\f$ or \f$X X^T\f$ where \f$X\f$ is
  * the input data matrix
  * \param len number of elements of input vector
  * \param r rank
  * \param n number of samples
  * \param nlambda
  * \param minl
  * \return vector of values for the regularization parameter
  */
template<typename T>
T* lambdaguesses(const T* eigvals, const int len, const int r, const int n, const int nlambda, const T minl)
{
    T* guesses = new T[nlambda];

    T* tmp = new T[len];
    copy(tmp, eigvals, len);

    std::sort(tmp, tmp+len);

    /*const*/ T lmin = tmp[len-r];
    const T lmax = tmp[len-1];

    delete[] tmp;

    T thr1 = std::min(lmin, minl*lmax);
    T thr2 = 200*static_cast<T>(sqrt(std::numeric_limits<T>::epsilon()));

    lmin = std::max(thr1, thr2);

    const T base = (lmax/lmin);
    const T den = nlambda-static_cast<T>(1.0);
    const T nT = (T)n;

    for(int i=0; i< nlambda; ++i)
        guesses[i] = (lmin * pow( base, ((T)i)/den ) )/nT;

    return guesses;
}

//template<typename T>
//T* mldivide(const T* A, const T* B,
//            const int a_rows, const int a_cols,
//            const int b_rows, const int b_cols,
//            int& res_rows, int& res_cols)
//{
//    // A\B = pinv(A)*B

//    int aInv_rows, aInv_cols;
//    T* A_inv = pinv(A, a_rows, a_cols, aInv_rows, aInv_cols);

//    T* ret = new T[aInv_rows*b_cols];
//    dot(A_inv, B, ret, aInv_rows, aInv_cols, b_rows, b_cols, aInv_rows, b_cols, CblasNoTrans, CblasNoTrans, CblasColMajor);

//    delete[] A_inv;

//    res_rows = aInv_rows;
//    res_cols = b_cols;

//    return ret;

//}

/**
  * Performs left division of squared matrices
  *
  * \param A first input matrix
  * \param B second input matrix. On exit B contains the result matrix
  * \param a_rows number of rows of matrix A
  * \param a_cols number of columns of matrix A
  * \param b_rows number of rows of matrix B
  * \param b_cols number of columns of matrix B
  * \param transA specifies whether to transpose A or not
  */
template<typename T>
void mldivide_squared(const T* A, T* B,
            const int a_rows, const int a_cols,
            const int b_rows, const int b_cols,
            const CBLAS_TRANSPOSE transA)
{
    if (a_cols != a_rows)
        throw gException("The input matrix A must be squared");

    const int lda = a_rows;
    const int ldb = b_rows;

    trsm(CblasLeft, CblasUpper, transA, CblasNonUnit, b_rows, b_cols, (T)1.0, A, lda, B, ldb);
}


/**
  * Template function to call BLAS *TRSM routines
  */
template <typename T>
void trsm(const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N, const T alpha, const T *A, const int lda, T *B, const int ldb);

/**
  * Computes singular value decomposition of an input matrix A such that  A = U*diag(S)*Vt.
  *
  * \param A Input matrix to be decomposed
  * \param U Matrix of the left singular vectors
  * \param S Vector containing the singular values of the decomposition
  * \param Vt transposed matrix of the right singular vectors
  * \param A_rows number of rows of matrix A
  * \param A_cols number of columns of matrix A
  * \param U_rows number of rows of matrix U
  * \param U_cols number of columns of matrix U
  * \param S_len number of elements of vector S
  * \param Vt_rows number of rows of matrix Vt
  * \param Vt_cols number of columns of matrix Vt
  * \param econ if true computes the "economy size" decomposition.
    If A_rows >= A_cols, then svd computes only the first A_cols columns of U and S length is A_cols.
    For A_rows < A_cols, only the first A_rows rows of Vt are computed and S length is A_rows.
  */
template <typename T>
void svd(const T* A, T*& U, T*& S, T*& Vt,
         const int A_rows, const int A_cols,
         int& U_rows, int& U_cols,
         int& S_len,
         int& Vt_rows, int& Vt_cols, bool econ = false) throw(gException)
{
    // A = U*S*Vt

    char jobu, jobvt;

    int m = A_rows;
    int n = A_cols;
    int k = std::min<int>(m, n);
    int ldvt;

    S = new T[k];
    S_len = k;

    U_rows = m;
    Vt_cols = n;

    if(econ)
    {
//        jobu = jobvt = 'S';
        if(m>n)
        {
            jobu = 'S';
            jobvt = 'A';
        }
        else if(m<n)
        {
            jobu = 'A';
            jobvt = 'S';
        }
        else
            jobu = jobvt = 'S';

        U = new T[m*k];
        Vt = new T[k*n];

        U_cols = k;
        Vt_rows = k;

        ldvt = k;
    }
    else
    {
        jobu = jobvt = 'A';

        U = new T[m*m];
        Vt = new T[n*n];

        U_cols = m;
        Vt_rows = n;

        ldvt = n;
    }

    //MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))

    int lda = A_rows;
    int ldu = m;
    int info, lwork = std::max<int>(3*k+std::max<int>(m,n), 5*k);
    T* work = new T[lwork];
    T* cpy = new T[m*n];
    copy(cpy, A, A_rows*A_cols);

    gesvd_(&jobu, &jobvt, &m, &n, cpy, &lda, S, U, &ldu, Vt, &ldvt, work, &lwork, &info);

    delete[] work;
    delete[] cpy;

    if(info != 0)
    {
        std::stringstream str;
        str << "SVD failed, error code " << info << std::endl;
        throw gException(str.str());
    }
}

/**
  * Template function to call BLAS *GESVD routines
  */
template<typename T>
int gesvd_(char *jobu, char *jobvt, int *m, int *n, T *a, int *lda, T *s, T *u, int *ldu, T *vt, int *ldvt, T *work, int *lwork, int *info);

/**
  * Generates a vector containing a random permutation of the values from start to start+n inclusive
  */
template<typename T>
void randperm(const unsigned long n, T* seq, bool generate = true, unsigned long start = 1)
{
    if(generate)
    {
        unsigned long val = start;
        for(T *it = seq, *end= seq+n; it != end; ++it)
            *it = val++;
    }

    for(unsigned long i=0; i<n; ++i)
        std::swap(seq[rand()%n], seq[rand()%n]);
}

/**
  * Generates a vector containing a copy of a row of an input matrix
  *
  * \param M input matrix
  * \param rows number of rows of the input matrix
  * \param cols number of columns of the input matrix
  * \param row_index index of the row to be extracted
  * \param row vector containing the extracted row. Length must be equal to \c cols
  */
template<typename T>
void getRow(const T* M, const int rows, const int cols, const int row_index , T* row)
{
    copy(row, M+row_index, cols, 1, rows);
}

/**
  * Template function to call BLAS *NRM2 routines
  */
template<typename T>
T nrm2(const int N, const T* X, const int incX);

/**
  * Template function to call BLAS *SCAL routines
  */
template<typename T>
void scal(const int N, const T alpha, T *X, const int incX);


/**
  * Template function to call BLAS *GEMV routines
  */
template<typename T>
void gemv(const CBLAS_TRANSPOSE TransA,
          const int M, const int N, const T alpha, const T *A, const int lda,
          const T *X, const int incX,
          const T beta, T *Y, const int incY);

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
  * Computes the eigenvalues/eigenvectors of a squared and symmetric input matrix.
  *
  * \param A input matrix. On exit it contains the orthonormal eigenvectors of the matrix A
  * \param L vector of eigenvalues in ascending order
  * \param A_rows_cols number of rows/columns of matrix A
  */
template<typename T>
void eig_sm(T* A, T* L, int A_rows_cols) throw (gException)
{
    char jobz = 'V';
    char uplo = 'L';
    int n = A_rows_cols, lda = n;
    int info, lwork = 4*n;
    T* work = new T[lwork];

    syev(&jobz, &uplo, &n, A, &lda, L, work, &lwork, &info );

    delete[] work;

    if(info != 0)
    {
        std::stringstream str;
        str << "Eigenvalues/eigenVectors computation failed, error code " << info << ";" << std::endl;
        throw gException(str.str());
    }
}

/**
  * In place computation of the exponential for each element of a vector
  *
  * \param v vector
  * \param length vector size
  */
template<typename T>
void exp(T* v, const int length)
{
    for(T *it = v, *end = v+length; it != end; ++it)
        *it = (T) std::exp(*it);
}

/**
  * Computes Euclidean distance between two vectors
  *
  * \param A first vector
  * \param B first vector
  * \param len number of element of vectors A and B
  * \param work work buffer of length >= \c len
  * \return Euclidean distance
  */
template<typename T>
T eucl_dist(const T* A, const T* B, const int len, T* work)
{
    copy(work, A, len);
    axpy(len, (T)-1.0, B, 1, work, 1);
    return nrm2(len, work, 1);
}

/**
  * Computes a vector D containing the Euclidean distances
  * between each pair of observations in the N-by-P data matrix A. Rows of
  * A correspond to observations, columns correspond to variables. D is a
  * 1-by-(N*(N-1)/2) vector, corresponding to the N*(N-1)/2 pairs of
  * observations in A.
  */
template<typename T>
void pdist(const T* A, const int N/*A_rows*/, const int P/*A_cols*/, T* D)
{
    T* work = new T[P];
    T* rowN = new T[P];
    T* rowNPlusOne = new T[P];
    T* D_it = D;

    for(int i=0; i<N-1; ++i)
    {
        getRow(A, N, P, i, rowN);

        for(int j=i+1; j<N; ++j)
        {
            getRow(A, N, P, j, rowNPlusOne);

            *D_it = eucl_dist(rowN,rowNPlusOne,P, work);
            ++D_it;
        }
    }

    delete [] work;
    delete [] rowN;
    delete [] rowNPlusOne;
}


/**
  * Reformats a distance matrix between upper triangular and square form.
  * If A is a vector as created by the pdist function,
  * converts A into a symmetric, square format, so that D(i,j) denotes the
  * distance between the i and j objects in the original data.
  * If A is a symmetric, square matrix with zeros along
  * the diagonal, creates a vector D containing the A elements below the
  * diagonal. D has the same format as the output from the PDIST function.
  */
template<typename T>
void squareform(const T* A, const int N/*A_rows*/, const int P/*A_cols*/, T* D, const int d_cols)
{
    if(d_cols!=1)
    {
        T* work = new T[P];
        T* rowN = new T[P];
        T* rowNPlusOne = new T[P];

        set(D, (T)0.0, N, N+1); // zeroes the diagonal

        for(int i=0; i<N; ++i)
        {
            copy(D+(i*N), D+i, i, 1, N); // copy from the other side

            if(i+1 < N)
            {
                getRow(A, N, P, i, rowN);

                for(int j=i+1; j<N; ++j)
                {
                    getRow(A, N, P, j, rowNPlusOne);
                    D[j + i*N] = eucl_dist(rowN, rowNPlusOne, P, work);
                }
            }
        }

        delete [] work;
        delete [] rowN;
        delete [] rowNPlusOne;

    }
    else
    {
        pdist(A, N, P, D);
    }
}

/**
  * Computes a vector containing the indices of all elements of an input vector that are equals to a given value
  *
  * \param V input vector
  * \param len number of elements of the input vector
  * \param value value for comparison
  * \param ind vector of the indices. Length must be >= len
  * \param ind_length total number of indices written into \c ind
  */
template<typename T>
void indicesOfEqualsTo(const T* V, const int len, const T value, unsigned long* ind, int& ind_length)
{
    unsigned long* r_it = ind;
    for(const T *it = V, *end = V+len; it != end; ++it)
    {
        if(eq(*it, value))
        {
            *r_it = it-V;
            ++r_it;
        }
    }

    ind_length = r_it - ind;
}

/**
  * Rounds an input value to the nearest integer
  */
template<typename T>
int round(const T value)
{
    if(eq(value, (T)0.0))
        return 0;

    return gt(value,(T)0.0)? ((int)(value+0.5)): ((int)(value-0.5));
}

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

/**
  * Computes an economy-size QR decomposition of an input matrix A so that A(:,E) = Q*R.
  * If m > n, only the first n columns of Q and the first n rows of R are computed.
  *
  * \param A input matrix
  * \param m number of rows of the input matrix
  * \param n number of columns of the input matrix
  * \param Q  Q matrix
  * \param R  upper triangular R matrix
  * \param E permutations vector
  */
template<typename T>
void qr_econ(const T* A, int m, int n, T* Q, T* R, int* E)
{
    // Q: mxmin(m,n)
    // R: min(m,n)xn
    // E: n

    T* Q_tmp = new T[m*n];
    copy(Q_tmp, A, m*n);

    int k = std::min(m, n);

    int lda = std::max(1,m);
    int* jpvt = E;
    set(jpvt, 0, n);

    T* tau = new T[k];

    T* work;
    int lwork = -1;
    int info;

    // query
    T qwork;
    geqp3( &m, &n, Q_tmp, &lda, jpvt, tau, &qwork, &lwork, &info);

    // exec
    lwork = static_cast<int>(qwork);
    work = new T[lwork];
    geqp3( &m, &n, Q_tmp, &lda, jpvt, tau, work, &lwork, &info);

    if(info != 0)
    {
        delete[] tau;
        delete[] work;
        delete[] Q_tmp;

        std::stringstream str;
        str << "QR factorization failed, error code " << info << std::endl;
        throw gException(str.str());
    }

    if(R != NULL)
    {
        for(int i=0; i<n; ++i)
        {
            copy(R + k*i, Q_tmp + m*i, i+1);

            set(R + (k*i) + i+1, (T)0.0, k - (i+1));
        }
    }

    // query
    lwork = -1;
    orgqr(&m, &k, &k, Q_tmp, &lda, tau, &qwork, &lwork, &info);

    if(info != 0)
    {
        delete[] tau;
        delete[] work;
        delete[] Q_tmp;

        std::stringstream str;
        str << "QR factorization failed, error code " << info << std::endl;
        throw gException(str.str());
    }

    //exec
    lwork = static_cast<int>(qwork);
    delete[] work;
    work = new T[lwork];
    orgqr(&m, &k, &k, Q_tmp, &lda, tau, work, &lwork, &info);

    copy(Q, Q_tmp, m*k);

    delete[] tau;
    delete[] work;
    delete[] Q_tmp;

}

/**
  * Generates a row vector of n points linearly spaced between and including a and b
  */
template<typename T>
void linspace(T a, T b, unsigned long n, T* res)
{
    unsigned long i = 0;

    if(n > 0)
    {
        const T coeff = (b - a)/(static_cast<T>(n)-1);
        for(i = 0; i < n-1; ++i)
            res[i] = a + (i*coeff);
    }

    res[i] = b;
}


/**
  * Generates a row vector containing the standard deviation of the elements of each column of an input matrix
  *
  * \param X input matrix
  * \param rows number of rows of the input matrix
  * \param cols number of columns of the input matrix
  * \param res output standard deviations vector
  * \param work work buffer of size >= cols+rows
  */
template<typename T>
void stdDev(const T* X, const int rows, const int cols, T* res, T* work)
{
    T* meanX = work;
    mean(X, meanX, rows, cols, cols);

    T* stdX = res;
    T* column = work+cols;

    for(int i=0; i< cols; ++i)
    {
        copy(column, X+(rows*i), rows);

        axpy(rows, (T)-1.0, meanX+i, 0, column, 1);

        stdX[i] = sqrt( pow(nrm2(rows, column, 1), 2) / (rows-1));
    }
}

/**
  * Returns the median value of a vector. This function does not preserve input vector.
  *
  * \param v input vector
  * \param length number of elements of the input vector
  */
template<typename T>
T median(T* v, const int length)
{
    std::sort(v, v+length);

    if(length%2)
        return v[length/2];
    else
        return static_cast<T>( (v[length/2] + v[(length/2)-1]) /2.0 );
}

/**
  * Computes the median values of the elements along different dimensions of a matrix
  *
  * \param M input matrix
  * \param rows number of rows of the input matrix
  * \param cols number of columns of the input matrix
  * \param res output median values vector, size == cols if dimension == 1  or size == rows if dimension == 2
  * \param work work buffer of size >= rows if dimension == 1  or >= cols if dimension == 2
  */
template<typename T>
void median(const T* M, const int rows, const int cols, const int dimension, T* res, T* work)
{
    switch(dimension)
    {
    case 1:
        for(int i=0; i<cols; ++i)
        {
            copy(work, M+(i*rows), rows);
            res[i] = median(work, rows);
        }
        break;
    case 2:
        for(int i=0; i<rows; ++i)
        {
            copy(work, M+i, cols, 1, rows);
            res[i] = median(work, cols);
        }
        break;
    default:
        throw gException(Exception_Illegal_Argument_Value);
    }
}

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

}

#endif // _GURLS_GMATH_H_
