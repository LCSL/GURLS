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
#include <cfloat>

#include "exports.h"
#include "gmat2d.h"
#include "gvec.h"
#include "exceptions.h"

#ifdef _WIN32
#pragma warning(disable : 4290)
#endif

#ifdef _ACML
#include <acml.h>
enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114};
enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142};
#elif defined _MKL
#include <mkl_cblas.h>
#else
extern "C" {
#include <cblas.h>
}
#endif

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
// ================= TODO: TO BE DISCUSSED =======================
// =========================================================
/* Replicate a vector n times along the columns (or along the rows if transpose==true)
   If x is a vector of length N, then the output is an N-by-N matrix whose columns (or rows) are copies of x.
*/
template <typename T>
gMat2D<T> repmat(const gVec<T>& x, unsigned long n, bool transpose = false){

    T dataA[n*x.size()];
    const T* datax = x.getData();
    unsigned long rows, cols;

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


template<typename T>
void set(T* buffer, const T value, const int size, const int incr);

template<typename T>
void set(T* buffer, const T value, const int size)
{
    for(T *it = buffer, *end = buffer+size; it != end; ++it)
        *it = value;
}

template<>
void GURLS_EXPORT set(float* buffer, const float value, const int size);

template<>
void GURLS_EXPORT set(double* buffer, const double value, const int size);


template<typename T>
void copy(T* dst, const T* src, const int size)
{
    memcpy(dst, src, size*sizeof(T));
}

template<>
GURLS_EXPORT void copy(float* dst, const float* src, const int size);

template<>
GURLS_EXPORT void copy(double* dst, const double* src, const int size);

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

template<>
GURLS_EXPORT void copy(float* dst, const float* src, const int size, const int dstIncr, const int srcIncr);

template<>
GURLS_EXPORT void copy(double* dst, const double* src, const int size, const int dstIncr, const int srcIncr);

//TODO optimize
template<typename T>
void copy_submatrix(T* dst, const T* src, const int src_Rows, const int sizeRows, const int sizeCols, unsigned long *indices_rows, unsigned long *indices_cols)
{
  int t=0;
  for (int j = 0; j < sizeCols; j++)
    for (int i = 0; i < sizeRows; i++)
    {
      dst[t] = src[ indices_rows[i] + indices_cols[j]*src_Rows ] ;
      ++t;
    }
}

template<typename T>
void subMatrixFromRows(const T* matrix, const int m_Rows, const int m_Cols, const unsigned long* rowsIndices, const int nIndices, T* submat, T* work)
{
    for(const unsigned long *it = rowsIndices, *end = rowsIndices+nIndices; it != end; ++it)
    {
        getRow(matrix, m_Rows, m_Cols, *it, work);

        copy(submat + (it-rowsIndices), work, m_Cols, nIndices, 1);
    }
}

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

template <typename T>
T* pinv(const T* A, const int rows, const int cols, int& res_rows, int& res_cols, T* RCOND = NULL);

template <typename T>
T* transpose_rm(const T* matrix, int rows, int cols)
{
    T* ret = new T[rows*cols];
    T* d1 = ret;
    const T* d0 = matrix;
    int N = cols;

    for (int c = 0; c < cols; ++c)
        for (int r = 0; r < rows; ++r)
            *d1++=*(d0+r*N+c);

    return ret;
}

template <typename T>
T* transpose_cm(const T* matrix, int rows, int cols)
{
    T* ret = new T[rows*cols];
    T* d1 = ret;
    const T* d0 = matrix;
    int N = rows;

    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            *d1++=*(d0+c*N+r);

    return ret;
}

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

template <typename T>
void dot(const T* A, const T* B, T* C,
         int A_rows, int A_cols,
         int B_rows, int B_cols,
         int C_rows, int C_cols,
         const /*enum*/ CBLAS_TRANSPOSE TransA,
         const /*enum*/ CBLAS_TRANSPOSE TransB,
         const /*enum*/ CBLAS_ORDER Order)
{

    bool transposeA = (TransA == CblasTrans || TransA == CblasConjTrans);
    bool transposeB = (TransB == CblasTrans || TransB == CblasConjTrans);

    if(transposeA)
        std::swap(A_rows, A_cols);

    if(transposeB)
        std::swap(B_rows, B_cols);

    if ((C_rows != A_rows) || (C_cols != B_cols))
        throw gException(gurls::Exception_Inconsistent_Size);

    const T alpha = 1.0;
    const T beta = 0.0;

    int lda, ldb, ldc;
    switch(Order)
    {
    case CblasColMajor:
        lda = transposeA? A_cols: A_rows;
        ldb = transposeB? B_cols: B_rows;
        ldc = C_rows;
        break;
    case CblasRowMajor:
        lda = transposeA? A_rows: A_cols;
        ldb = transposeB? B_rows: B_cols;
        ldc = C_cols;
        break;
    default:
        lda = ldb = ldc = 0;
        assert(0);
    }

    // C = alpha*A*B + beta*C
    gemm(Order, TransA, TransB, C_rows, C_cols, A_cols, alpha, A, lda, B, ldb, beta, C, ldc);

}

template<typename T>
T* cholesky(const T* matrix, const int rows, const int cols, bool upper = true)
{
    T* ret = new T[rows*cols];
    copy(ret, matrix, rows*cols);

    int LDA = rows;
    int nc = cols;
    char UPLO = upper? 'U':'L';
    int info;

    potrf_(&UPLO, &nc, ret, &LDA, &info);

    if(info != 0)
    {
        std::stringstream str;
        str << "Cholesky factorization failed, error code " << info << ";" << std::endl;
        throw gException(str.str());
    }

    clearLowerTriangular(ret, rows, cols);

    return ret;
}

template<typename T>
int potrf_(char *UPLO, int *n, T *a, int *lda , int *info);

template <typename T>
void eig(const T* A, T* &V, T* &Wr, T* &Wi, int A_rows, int A_cols) throw (gException);

template <typename T>
void eig(const T* A, T* &V, T* &W, int A_rows, int A_cols) throw (gException)
{
    T* tmp;
    try
    {
        eig(A, V, W, tmp, A_rows, A_cols);
    }
    catch(gException &ex)
    {
        delete[] tmp;
        throw(ex);
    }
    delete[] tmp;
}

template <typename T>
void setReciprocal(T* matrix, const int len)
{
    const T one = static_cast<T>(1.0);

    for(T *it = matrix, *end = matrix+len; it != end; ++it)
        *it = one / *it;
}

template <typename T>
T* diag(T* vector, const int len)
{
    T* ret = new T[len*len];

    const T zero = static_cast<T>(0.0);

    set(ret, zero, len*len);
    copy(ret, vector, len, len+1, 1);

    return ret;
}


template <typename T>
T mul1(T a, T b)
{
    return a*b;
}

template <typename T>
T div(T a, T b)
{
    return a/b;
}

template <typename T>
T diff(T a, T b)
{
    return a-b;
}

template<typename T>
bool eq(T val1, T val2)
{
    return val1 == val2;
}

template<typename T>
bool gt(T a, T b);

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

template <typename T>
void mult(const T* A, const T* B, T* result, const int len)
{
    binOperation<T>(A, B, result, len, &mul1);
}

template <typename T>
void rdivide(const T* A, const T* B, T* result, const int len)
{
    binOperation<T>(A, B, result, len, &div);
}

template <typename T>
void sum(const T* A, T* result, const int A_rows, const int A_cols, const int res_length) throw (gException)
{
    if(A_cols != res_length)
        throw gException("Sum: vector lengths mismatch");

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

template <typename T>
T sumv(const T* V, const int len, T* work) throw (gException)
{
    set(work, (T)1.0, len);

    gemv(CblasColMajor, CblasNoTrans, 1, len, (T)1.0, work, 1, V, 1, (T)0.0, work+len, 1);

    return work[len];
}

template <typename T>
void mean(const T* A, T* result, const int A_rows, const int A_cols, const int res_length) throw (gException)
{
    sum(A, result, A_rows, A_cols, res_length);

    for(T *it = result, *end = result+A_cols; it != end; ++it)
        *it /= A_rows;
}

//min along rows
template <typename T>
void argmin(const T* A, unsigned long* result, const int A_rows, const int A_cols, const int res_length) throw (gException)
{
    if(A_cols != res_length)
        throw gException("Sum: vector lengths mismatch");

    const T *a_it = A;

    for(unsigned long *r_it = result, *r_end = result+A_cols; r_it != r_end; ++r_it, a_it += A_rows)
        *r_it = (std::min_element(a_it, a_it+A_rows) - a_it);
}

template <typename T>
T* copyLocations(const unsigned long* locs, const T* src, const int locs_len, const int src_len)
{
    T* v = new T[locs_len];
    T* ptr_v = v;

    int val;
    for(const unsigned long* l_it = locs, *l_end=locs+locs_len; l_it != l_end; ++l_it, ++ptr_v)
    {
        val = *l_it;
        if((val < 0) || (val > src_len))
            throw gException(gurls::Exception_Index_Out_of_Bound);

        *ptr_v = src[val];
    }

    return v;
}

template <typename T>
void axpy(const int N, const T alpha, const T *X, const int incX, T *Y, const int incY);

//template <typename T>
//void axpby(const int N, const T alpha, const T *X, const int incX, const T beta, T *Y, const int incY);

template <typename T>
T dot(const int N, const T *X, const int incX, const T *Y, const int incY);

template<typename T>
void rls_eigen(const T* Q, const T* L, const T* Qty, T* C, const T lambda, const int n,
             const int Q_rows, const int Q_cols,
             const int L_length,
             const int Qty_rows, const int Qty_cols)//  throw (gException)
{
    //function C = rls_eigen(Q,L,QtY,lambda,n)

    //sQ = size(Q,1); -> Q_rows

    //L = L + n*lambda;
    T* L1 = new T[L_length];

    set(L1, n*lambda , L_length);

    axpy(L_length, (T)1.0, L, 1, L1, 1);

    //L = L.^(-1);
    setReciprocal(L1, L_length);

    //L = spdiags(L,0,sQ,sQ);
    T* QL = new T[Q_rows*L_length];

    copy(QL, Q, Q_rows*Q_cols);
    for(int i=0; i< Q_cols; ++i)
        scal(Q_rows, L1[i], QL+(Q_rows*i), 1);


    //C = (Q*L)*QtY;
//    T* C = new T[Q_rows*Qty_cols];
    dot(QL, Qty, C, Q_rows, L_length, Qty_rows, Qty_cols, Q_rows, Qty_cols, CblasNoTrans, CblasNoTrans, CblasColMajor);

    delete[] L1;
//    delete[] Ldiag;
    delete[] QL;

//    return C;
}


template<typename T>
void GInverseDiagonal(const T* Q, const T* L, const T* lambda, T* Z,
                    const int Q_rows, const int Q_cols,
                    const int L_length, const int lambda_length)
{
    //function Z = GInverseDiagonal( Q, L, lambda )

    //n = size(Q, 1); -> Q_rows
    //t = size(lambda, 2); -> lambda_length

    //Z = zeros(n, t);

    //D = Q.^(2);
    const int Q_size = Q_rows*Q_cols;
    T* D = new T[Q_size];
    //copy(D, Q, Q_size, 1, 1);
    mult(Q, Q, D, Q_size);

    T* d = new T[L_length];
//    T* Dd = new T[Q_rows /* *1 */];

    //for i = 1 : t
    for(int i=0; i<lambda_length; ++i)
    {
//    d = L + (n*lambda(i));
        set(d, Q_rows*lambda[i] , L_length);
        axpy(L_length, (T)1.0, L, 1, d, 1);

//    d  = d.^(-1);
        setReciprocal(d, L_length);

//    Z(:,i) = D*d;
        gemv(CblasColMajor, CblasNoTrans, Q_rows, Q_cols, (T)1.0, D, Q_cols, d, 1, (T)0.0, Z+(i*lambda_length), 1);
    }

    delete[] d;
    delete[] D;
}


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

template<typename T>
T* compare(const T* vector1, const T* vector2, const int size, bool(*pred)(T,T))
{
    const T* v1_end = vector1+size;

    T* ret = new T[size];
    T* r_it = ret;

    for(const T *v1_it = vector1, *v2_it = vector2; v1_it != v1_end; ++v1_it, ++v2_it, ++r_it)
    {
        *r_it = static_cast<T>((*pred)(*v1_it, *v2_it)? 1.0: 0.0);
    }

    return ret;
}

template<typename T>
T* compare(const T* vector1, const T val, const int size, bool(*pred)(T,T))
{
    const T* v1_end = vector1+size;

    T* ret = new T[size];
    T* r_it = ret;

    for(const T *v1_it = vector1; v1_it != v1_end; ++v1_it, ++r_it)
    {
        *r_it = (*pred)(*v1_it, val)? (T)1.0: (T)0.0;
    }

    return ret;
}

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
        throw gException(gurls::Exception_Illegat_Argument_Value);
    }

    for(unsigned long *r_it = ind, *r_end = ind+m_cols; r_it != r_end; ++r_it, m_it += m_rows)
        *r_it = (std::max_element(m_it, m_it+m_rows) - m_it);

}

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
    T thr2 = 200*static_cast<T>(sqrt(DBL_EPSILON));

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

    trsm(CblasColMajor, CblasLeft, CblasUpper, transA, CblasNonUnit, b_rows, b_cols, (T)1.0, A, lda, B, ldb);
}

template <typename T>
void trsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N, const T alpha, const T *A, const int lda, T *B, const int ldb);

template <typename T>
void svd(const T* A, T*& U, T*& W, T*& Vt,
         const int A_rows, const int A_cols,
         int& U_rows, int& U_cols,
         int& W_len,
         int& Vt_rows, int& Vt_cols) throw(gException)
{
    // A = U*S*Vt

    char jobu = 'S', jobvt = 'S';
    int m = A_rows;
    int n = A_cols;
    int k = std::min<int>(m, n);

    W = new T[k];
    U = new T[m*k];
    Vt = new T[k*n];

    U_rows = m;
    U_cols = k;

    W_len = k;

    Vt_rows = k;
    Vt_cols = n;
//MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
    int lda = A_cols;
    int ldu = k;
    int ldvt = n;
    int info, lwork = std::max<int>(3*k+std::max<int>(m,n), 5*k);
    T* work = new T[lwork];
    T* cpy = new T[m*n];
    copy(cpy, A, A_rows*A_cols);

    gesvd_(&jobu, &jobvt, &n, &m, cpy, &lda, W, U, &ldu, Vt, &ldvt, work, &lwork, &info);

    delete[] work;
    delete[] cpy;

    if(info != 0)
    {
        std::stringstream str;
        str << "SVD failed, error code " << info << std::endl;
        throw gException(str.str());
    }
}

template<typename T>
int gesvd_(char *jobu, char *jobvt, int *m, int *n, T *a, int *lda, T *s, T *u, int *ldu, T *vt, int *ldvt, T *work, int *lwork, int *info);

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

template<typename T>
void getRow(const T* M, const int rows, const int cols, const int row_index , T* row)
{
    copy(row, M+row_index, cols, 1, rows);
}

template<typename T>
T nrm2(const int N, const T* X, const int incX);

template<typename T>
void scal(const int N, const T alpha, T *X, const int incX);

template<typename T>
void gemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
          const int M, const int N, const T alpha, const T *A, const int lda,
          const T *X, const int incX,
          const T beta, T *Y, const int incY);

template<typename T>
void syev( char* jobz, char* uplo, int* n, T* a, int* lda, T* w, T* work, int* lwork, int* info);

template<typename T>
void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
          const int M, const int N, const int K, const T alpha, const T *A, const int lda,
          const T *B, const int ldb,
          const T beta, T *C, const int ldc);

template<typename T>
void eig_sm(T* A, T* L, int A_rows, int A_cols) throw (gException)
{
    if (A_cols != A_rows)
        throw gException("The input matrix A must be squared");

    char jobz = 'V';
    char uplo = 'L';
    int n = A_cols, lda = A_rows;
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

template<typename T>
void exp(T* v, const int length)
{
    for(T *it = v, *end = v+length; it != end; ++it)
        *it = (T) std::exp(*it);
}

template<typename T>
T eucl_dist(const T* A, const T* B, const int len, T* work)
{
    copy(work, A, len);
    axpy(len, (T)-1.0, B, 1, work, 1);
    return nrm2(len, work, 1);

}

// Y = PDIST(X) returns a vector Y containing the Euclidean distances
// between each pair of observations in the N-by-P data matrix X.  Rows of
// X correspond to observations, columns correspond to variables.  Y is a
// 1-by-(N*(N-1)/2) row vector, corresponding to the N*(N-1)/2 pairs of
// observations in X.
// D allocated
template<typename T>
void pdist(const T* A, const int N/*A_rows*/, const int P/*A_cols*/, T* D)
{
//    int d_len = N*(N-1)/2;
    T* work = new T[P];
    T* rowN = new T[P];
    T* rowNPlusOne = new T[P];
    //D = new T[d_len];
    int cont = 0;
    for(int i=0;i<N;i++)
    {
        copy< T >(rowN,A + i,P,1,N);
        //    getRow< T >(A,N,P,i,rowN);
        for(int j=i+1;j<N;j++)
        {
            copy< T >(rowNPlusOne,A + j,P,1,N);
            //      getRow< T >(A,N,P,j,rowNPlusOne);
            D[cont] = eucl_dist<T>(rowN,rowNPlusOne,P, work);
            cont++;
        }
    }
    delete [] work;
    delete [] rowN;
    delete [] rowNPlusOne;
}

//  SQUAREFORM Reformat a distance matrix between upper triangular and square form.
//     Z = SQUAREFORM(Y), if Y is a vector as created by the PDIST function,
//     converts Y into a symmetric, square format, so that Z(i,j) denotes the
//     distance between the i and j objects in the original data.
//   Y = SQUAREFORM(Z), if Z is a symmetric, square matrix with zeros along
//     the diagonal, creates a vector Y containing the Z elements below the
//     diagonal.  Y has the same format as the output from the PDIST function.
template<typename T>
void squareform(const T* A, const int N/*A_rows*/, const int P/*A_cols*/, T* D, const int d_cols)
{
    if(d_cols!=1)
    {
        //D = new T[N*N];
        T* work = new T[P];
        T* rowN = new T[P];
        T* rowNPlusOne = new T[P];
        for(int i=0;i<N;i++)
        {
            copy< T >(rowN,A + i,P,1,N);
            for(int j=0;j<=N;j++)
            {
                if(i!=j)
                {
                    copy< T >(rowNPlusOne,A + j,P,1,N);
                    D[j + i*N] = eucl_dist< T >(rowN,rowNPlusOne,P, work);
                }
                else
                    D[j + i*N]=0;
            }
        }
        delete [] work;
        delete [] rowN;
        delete [] rowNPlusOne;
    }
    else
    {
        //   D = new T[N*(N-1)/2];
        pdist<T>(A, N, P, D);
    }
}

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

template<typename T>
int round(const T value)
{
    if(eq(value, (T)0.0))
        return 0;

    return gt(value,(T)0.0)? ((int)(value+0.5)): ((int)(value-0.5));
}

}
#endif // _GURLS_GMATH_H_
