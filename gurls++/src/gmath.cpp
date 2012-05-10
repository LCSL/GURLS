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


extern "C" {

int sgetrf_(int *m, int *n, float *a, int *lda, int *ipiv, int *info);
int sgesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info);
int sgeev_(char *jobvl, char *jobvr, int *n, float *a, int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);
int sgetri_(int *n, float *a, int* lda, int *ipiv, float*work, int *lwork, int *info);
int spotrf_(char *UPLO, int *n, float *a, int *lda , int *info);
int sgelss_( int *m, int *n, int* nrhs, float *a, int *lda, float* b, int *ldb, float *s, float *rcond, int *rank, float *work, int *lwork, int *info);
int ssyev_( char* jobz, char* uplo, int* n, float* a, int* lda, float* w, float* work, int* lwork, int* info );
int dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

}

namespace gurls {

template <>
float dot(const gVec<float>& x, const gVec<float>& y) {

    if ( x.getSize() != y.getSize() ) {
        throw gException(gurls::Exception_Inconsistent_Size);
    }

    return ( cblas_sdot (x.getSize(), x.getData(), 1, y.getData(), 1) );

}

template <>
double dot(const gVec<double>& x, const gVec<double>& y) {

    if ( x.getSize() != y.getSize() ) {
        throw gException(gurls::Exception_Inconsistent_Size);
    }

    return ( cblas_ddot (x.getSize(), x.getData(), 1, y.getData(), 1) );

}


// ============ OUT-OF-PLACE MATRIX MULTIPLICATION ==================

template <>
void dot(const gMat2D<float>& A, const gMat2D<float>& B, gMat2D<float>& C) {

    if ((C.rows() != A.rows()) || (C.cols() != B.cols())) {
        throw gException(gurls::Exception_Inconsistent_Size);
    }

    // C = alpha*A*B + beta*C
    const float alpha = 1.;
    const float beta = 0.;
    const /*enum*/ CBLAS_TRANSPOSE TransA = CblasNoTrans;
    const /*enum*/ CBLAS_TRANSPOSE TransB = CblasNoTrans;
    const /*enum*/ CBLAS_ORDER Order = CblasRowMajor;
    cblas_sgemm(Order, TransA, TransB, C.rows(), C.cols(), A.cols(), alpha,
                A.getData(), A.cols(), B.getData(), B.cols(), beta, C.getData(), C.cols());

}


template <>
void dot(const gMat2D<double>& A, const gMat2D<double>& B, gMat2D<double>& C) {

    if ((C.rows() != A.rows()) || (C.cols() != B.cols())) {
        throw gException("Inconsistent matrix dimensions.");
    }

    // C = alpha*A*B + beta*C
    const double alpha = 1.;
    const double beta = 0.;
    const /*enum*/ CBLAS_TRANSPOSE TransA = CblasNoTrans;
    const /*enum*/ CBLAS_TRANSPOSE TransB = CblasNoTrans;
    const /*enum*/ CBLAS_ORDER Order = CblasRowMajor;
    cblas_dgemm(Order, TransA, TransB, C.rows(), C.cols(), A.cols(), alpha,
                A.getData(), A.cols(), B.getData(), B.cols(), beta, C.getData(), C.cols());

}



// ============ OUT-OF-PLACE MATRIX-VECTOR MULTIPLICATION ==================

template <>
void dot(const gMat2D<float>& A, const gVec<float>& x, gVec<float>& y){

    if ( (A.cols() != x.getSize()) ||  (A.rows() != y.getSize())){
        throw gException("Inconsistent matrix dimensions.");
    }

    // y = alpha*A*x + beta*y
    const float alpha = 1.;
    const float beta = 0.;
    const /*enum*/ CBLAS_TRANSPOSE TransA = CblasNoTrans;
    const /*enum*/ CBLAS_ORDER Order = CblasRowMajor;

    cblas_sgemv(Order, TransA, A.rows(), A.cols(),alpha, A.getData(), A.cols(),
                x.getData(), 1, beta, y.getData(), 1);
}

template <>
void dot(const gMat2D<double>& A, const gVec<double>& x, gVec<double>& y){
    if ( (A.cols() != x.getSize()) ||  (A.rows() != y.getSize())){
        throw gException("Inconsistent matrix dimensions.");
    }

    // y = alpha*A*x + beta*y
    const double alpha = 1.;
    const double beta = 0.;
    const /*enum*/ CBLAS_TRANSPOSE TransA = CblasNoTrans;
    const /*enum*/ CBLAS_ORDER Order = CblasRowMajor;

    cblas_dgemv(Order, TransA, A.rows(), A.cols(),alpha, A.getData(), A.cols(),
                x.getData(), 1, beta, y.getData(), 1);

}

template <>
void lu(gMat2D<float>& A, gVec<int>& pv) {

    unsigned int k = std::min<unsigned int>(A.cols(), A.rows());
    if (pv.getSize() != k) {
        throw gException("The lenghth of pv must be equal to the minimun dimension of A");
    }
    int info;
    int m = A.rows();
    int n = A.cols();
    int lda = A.cols();
    sgetrf_(&m, &n, A.getData(), &lda, pv.getData(), &info);

}

template <>
void lu(gMat2D<float>& A) {

    gVec<int> pv(std::min<int>(A.cols(), A.rows()));
    lu(A, pv);
}

template <>
void inv(const gMat2D<float>& A, gMat2D<float>& Ainv, InversionAlgorithm alg){
    Ainv = A;
    int k = std::min<int>(Ainv.cols(), Ainv.rows());
    int info;
    int* ipiv = new int[k];
    int m = Ainv.rows();
    int n = Ainv.cols();
    int lda = Ainv.cols();
    float* work = new float[n];

    sgetrf_(&m, &n, Ainv.getData(), &lda, ipiv, &info);

    sgetri_(&m, Ainv.getData(), &lda, ipiv, work, &n, &info);
    delete[] ipiv;
    delete[] work;
}


template <>
void pinv(const gMat2D<float>& A, gMat2D<float>& Ainv, float RCOND){

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

    int M = A.rows();
    int N = A.cols();

    // The following step is required because we are currently storing
    // the matrices using a column-major order while LAPACK's
    // routines require row-major ordering
    float* a = new float[M*N];
    const float* ptr_A = A.getData();
    float* ptr_a = a;
    for (int j = 0; j < N ; j++){
        for (int i = 0; i < M ; i++){
            *ptr_a++ = *(ptr_A+i*N+j);
        }
    }
    int LDA = M;
    int LDB = std::max(M, N);
    int NRHS = LDB;

    float *b = new float[LDB*NRHS], *b_ptr = b;
    for (int i = 0; i < LDB*NRHS; i++){
        *b_ptr++=0.f;
    }
    b_ptr = b;
    for (int i = 0; i < std::min(LDB, NRHS); i++, b_ptr+=(NRHS+1)){
        *b_ptr = 1.f;
    }

    float* S = new float[std::min(M,N)];
    float condnum = 0.f; // The condition number of A in the 2-norm = S(1)/S(min(m,n)).

    if (RCOND < 0){
        RCOND = 0.f;
    }
    int RANK = -1; // std::min(M,N);
    int LWORK = -1; //2 * (3*LDB + std::max( 2*std::min(M,N), LDB));
    float* WORK = new float[1];
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
    int res = sgelss_( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);
    LWORK = WORK[0];
    delete [] WORK;
    WORK = new float[LWORK];

    res = sgelss_( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);
    // TODO: check INFO on exit
    condnum = S[0]/(S[std::min(M, N)]-1);


    gMat2D<float> *tmp = new gMat2D<float>(b, LDB, LDB, false);

    float *ainv = new float[N*M];
    float* ptr_b = ainv;
    float* ptr_B = b;
    for (int i = 0; i < N ; i++){
        for (int j = 0; j < M ; j++){
            *(ptr_b+i*M+j) = *(ptr_B+j*NRHS+i);
        }
    }
    Ainv = * new gMat2D<float>(ainv, N, M, true);

//	gMat2D<float> *tmp = new gMat2D<float>(b, LDB, NRHS, false);
//	gMat2D<float> *tmp1 = new gMat2D<float>(NRHS, LDB);
//	tmp->transpose(*tmp1);
//	Ainv = * new gMat2D<float>(tmp1->getData(), N, M, true);
//	std::cout << "A = " << std::endl << A << std::endl;
//	std::cout << "pinv(A) = " << std::endl << Ainv << std::endl;
    delete [] S;
    delete [] WORK;
    delete [] a;
//	delete tmp, tmp1;
    delete [] b;
}

template <>
void svd(const gMat2D<float>& A, gMat2D<float>& U, gVec<float>& W, gMat2D<float>& Vt) {

    char jobu = 'S', jobvt = 'S';
    int m = A.rows();
    int n = A.cols();
    int k = std::min<int>(m, n);

    if (W.getSize() < k) {
        throw gException("The length of vector W must be at least equal to the minimum dimension of the input matrix A");
    }
    if (U.rows() < m || U.cols() < k) {
        throw gException("Please check the dimensions of the matrix U where to store the singular vectors");
    }
    if (Vt.rows() < k || Vt.cols() < n) {
        throw gException("Please check the dimensions of the matrix Vt where to store the rigth singular vectors");
    }

    int lda = A.cols();
    int ldu = U.cols();
    int ldvt = Vt.cols();
    int info, lwork = std::max<int>(3*k+std::max<int>(m,n), 5*k);
    float* work = new float[lwork];
    float* copy = new float[m*n];
    A.asarray(copy, m*n);
    sgesvd_(&jobu, &jobvt, &n, &m, copy, &lda, W.getData(), Vt.getData(), &ldvt, U.getData(), &ldu, work, &lwork, &info);
    delete[] work;
    delete[] copy;
}


template <>
void eig(const gMat2D<float>& A, gMat2D<float>& V, gVec<float>& Wr, gVec<float>& Wi) {

    if (A.cols() != A.rows()) {
        throw gException("The input matrix A must be squared");
    }

    char jobvl = 'N', jobvr = 'V';
    int n = A.cols(), lda = A.cols(), ldvl = 1, ldvr = A.cols();
    int info, lwork = 4*n;
    float* work = new float[lwork];
    gMat2D<float> Atmp = A;
    gMat2D<float> Vtmp = V;
    sgeev_(&jobvl, &jobvr, &n, Atmp.getData(), &lda, Wr.getData(), Wi.getData(), NULL, &ldvl, Vtmp.getData(), &ldvr, work, &lwork, &info);
    Vtmp.transpose(V);
    delete[] work;
}


template <>
void eig(const gMat2D<float>& A, gMat2D<float>& V, gVec<float>& W){
    gVec<float> tmp(W.getSize());
    tmp = 0;
    eig(A, V, W, tmp);
}

template <>
void eig(const gMat2D<float>& A, gVec<float>& Wr, gVec<float>& Wi) {

    if (A.cols() != A.rows()) {
        throw gException("The input matrix A must be squared");
    }

    char jobvl = 'N', jobvr = 'N';
    int n = A.cols(), lda = A.cols(), ldvl = 1, ldvr = 1;
    int info, lwork = 4*n;
    float* work = new float[lwork];
    sgeev_(&jobvl, &jobvr, &n, const_cast<gMat2D<float>&>(A).getData(), &lda, Wr.getData(), Wi.getData(), NULL, &ldvl, NULL, &ldvr, work, &lwork, &info);
    delete[] work;
}

template <>
void eig(const gMat2D<float>& A, gVec<float>& W){
    gVec<float> tmp = W;
    eig(A, W, tmp);
}



template <>
void cholesky(const gMat2D<float>& A, gMat2D<float>& L, bool upper){

    typedef float T;
    L = A;

    int LDA = A.rows();
    int n = A.cols();
    char UPLO = upper? 'U' : 'L';
    int info;

    spotrf_(&UPLO,&n, L.getData(),&LDA,&info);

    // This is required because we adopted a column major order to store the
    // data into matrices
    gMat2D<T> tmp(L.rows(), L.cols());
    if (!upper){
        L.uppertriangular(tmp);
    } else {
        L.lowertriangular(tmp);
    }
    tmp.transpose(L);

}



template<>
void set(float* buffer, const float value, const int size, const int incr)
{
    cblas_scopy(size, &value, 0, buffer, incr);
}

template<>
void set(float* buffer, const float value, const int size)
{
    set(buffer, value, size, 1);
}

template<>
void set(double* buffer, const double value, const int size, const int incr)
{
    cblas_dcopy(size, &value, 0, buffer, incr);
}

template<>
void set(double* buffer, const double value, const int size)
{
    set(buffer, value, size, 1);
}


template<>
void copy(float* dst, const float* src, const int size, const int dstIncr, const int srcIncr)
{
    cblas_scopy(size, src, srcIncr, dst, dstIncr);
}

template<>
void copy(float* dst, const float* src, const int size)
{
    cblas_scopy(size, src, 1, dst, 1);
}

template<>
void copy(double* dst, const double* src, const int size, const int dstIncr, const int srcIncr)
{
    cblas_dcopy(size, src, srcIncr, dst, dstIncr);
}

template<>
void copy(double* dst, const double* src, const int size)
{
    cblas_dcopy(size, src, 1, dst, 1);
}

template<>
float* pinv(const float* A, int rows, int cols, int& res_rows, int& res_cols, float* RCOND)
{
    int M = rows;
    int N = cols;

    float* a = new float[rows*cols];
    copy(a, A, rows*cols);

    int LDA = M;
    int LDB = std::max(M, N);
    int NRHS = LDB;


    // float* b = eye(LDB).getData()

    const int b_size = LDB*NRHS;
    float *b = new float[LDB*NRHS];
    set(b, 0.f, b_size);
    set(b, 1.f, std::min(LDB, NRHS), NRHS+1);

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
    int res = sgelss_( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &rcond, &RANK, WORK, &LWORK, &INFO);
    LWORK = WORK[0];
    delete [] WORK;
    WORK = new float[LWORK];

    res = sgelss_( &M, &N, &NRHS, a, &LDA, b, &LDB, S, &rcond, &RANK, WORK, &LWORK, &INFO);

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

template<>
void dot(const float* A, const float* B, float* C,
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

    const float alpha = 1.0;
    const float beta = 0.0;

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
//    default:
//        assert(0);
    }

    // C = alpha*A*B + beta*C
    cblas_sgemm(Order, TransA, TransB, C_rows, C_cols, A_cols, alpha,
                A, lda, B, ldb, beta, C, ldc);

}

template<>
float* cholesky(const float* matrix, const int rows, const int cols, bool upper)
{
    float* ret = new float[rows*cols];
    copy(ret, matrix, rows*cols);

    int LDA = rows;
    int nc = cols;
    char UPLO = upper? 'U':'L';
    int info;

    spotrf_(&UPLO, &nc, ret, &LDA, &info);

    if(info != 0)
    {
        std::stringstream str;
        str << "Cholesky factorization failed, error code " << info << ";" << std::endl;
        throw gException(str.str());
    }

    clearLowerTriangular(ret, rows, cols);

    return ret;
}

template<>
void eig(const float* A, float* &V, float* &Wr, float* &Wi, int A_rows, int A_cols) throw (gException)
{
    if (A_cols != A_rows)
        throw gException("The input matrix A must be squared");

    const int size = A_rows*A_cols;
    float* Atmp = new float[size];
    copy(Atmp, A, size);

    V = new float[size];
    Wr = new float[A_rows];
    Wi = new float[A_rows];

    char jobvl = 'N', jobvr = 'V';
    int n = A_cols, lda = A_cols, ldvl = 1, ldvr = A_cols;
    int info, lwork = 4*n;
    float* work = new float[lwork];

    sgeev_(&jobvl, &jobvr, &n, Atmp, &lda, Wr, Wi, NULL, &ldvl, V, &ldvr, work, &lwork, &info);

    delete[] Atmp;
    delete[] work;

    if(info != 0)
    {
        std::stringstream str;
        str << "Eigenvalues/eigenVectors computation failed, error code " << info << ";" << std::endl;
        throw gException(str.str());
    }
}

template<>
void eig(const float* A, float* &V, float* &W, int A_rows, int A_cols) throw (gException)
{
    float* tmp;
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

template<>
void axpy(const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
{
    cblas_saxpy(N, alpha, X, incX, Y, incY);
}

template<>
void axpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY)
{
    cblas_daxpy(N, alpha, X, incX, Y, incY);
}

//template <>
//void axpby(const int N, const float alpha, const float *X, const int incX, const float beta, float *Y, const int incY)
//{
//    cblas_saxpby(N, alpha, X, incX, beta, Y, incY);
//}

//template <>
//void axpby(const int N, const double alpha, const double *X, const int incX, const double beta, double *Y, const int incY)
//{
//    cblas_daxpby(N, alpha, X, incX, beta, Y, incY);
//}

template <>
float dot(const int N, const float *X, const int incX, const float *Y, const int incY)
{
    cblas_sdot(N, X, incX, Y, incY);
}

template <>
double dot(const int N, const double *X, const int incX, const double *Y, const int incY)
{
    cblas_ddot(N, X, incX, Y, incY);
}

template<>
bool eq(double val1, double val2)
{
    return (val1 >= val2-DBL_EPSILON && val1 <= val2+DBL_EPSILON );
}

template<>
bool eq(float val1, float val2)
{
    return ( val1 >= val2-FLT_EPSILON && val1 <= val2+FLT_EPSILON );
}

template <>
void svd(const float* A, float*& U, float*& W, float*& Vt,
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

    W = new float[k];
    U = new float[m*k];
    Vt = new float[k*n];

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
    float* work = new float[lwork];
    float* cpy = new float[m*n];
    copy(cpy, A, A_rows*A_cols);

    sgesvd_(&jobu, &jobvt, &n, &m, cpy, &lda, W, U, &ldu, Vt, &ldvt, work, &lwork, &info);

    delete[] work;
    delete[] cpy;

    if(info != 0)
    {
        std::stringstream str;
        str << "SVD failed, error code " << info << std::endl;
        throw gException(str.str());
    }
}

void randperm(const unsigned int n, unsigned int* seq)
{
    unsigned int val = 0;
    for(unsigned int *it = seq, *end= seq+n; it != end; ++it)
        *it = ++val;

    for(int i=0; i<n; ++i)
        std::swap(seq[rand()%n], seq[rand()%n]);
}

template<>
float nrm2(const int N, const float* X, const int incX)
{
    cblas_snrm2(N, X, incX);
}

template<>
double nrm2(const int N, const double* X, const int incX)
{
    cblas_dnrm2(N, X, incX);
}

template<>
void scal(const int N, const float alpha, float *X, const int incX)
{
    cblas_sscal(N, alpha, X, incX);
}

template<>
void scal(const int N, const double alpha, double *X, const int incX)
{
    cblas_dscal(N, alpha, X, incX);
}

template<>
void gemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
          const int M, const int N, const float alpha, const float *A,
          const int lda, const float *X, const int incX,
          const float beta, float *Y, const int incY)
{
    cblas_sgemv(order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

template<>
void gemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
          const int M, const int N, const double alpha, const double *A,
          const int lda, const double *X, const int incX,
          const double beta, double *Y, const int incY)
{
    cblas_dgemv(order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}

template<>
void syev( char* jobz, char* uplo, int* n, float* a, int* lda, float* w, float* work, int* lwork, int* info)
{
    ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

template<>
void syev( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info)
{
    dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

template <>
void trsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N, const float alpha, const float *A, const int lda, float *B, const int ldb)
{
    cblas_strsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template <>
void trsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N, const double alpha, const double *A, const int lda, double *B, const int ldb)
{
    cblas_dtrsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

}
