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
#include <cblas.h>

int sgetrf_(int *m, int *n, float *a, int *lda, int *ipiv, int *info);
int sgesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info);
int sgeev_(char *jobvl, char *jobvr, int *n, float *a, int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);
int sgetri_(int *n, float *a, int* lda, int *ipiv, float*work, int *lwork, int *info);
int spotrf_(char *UPLO, int *n, float *a, int *lda , int *info);
int sgelss_( int *m, int *n, int* nrhs, float *a, int *lda, float* b, int *ldb, float *s, float *rcond, int *rank, float *work, int *lwork, int *info);

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
	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
	const enum CBLAS_ORDER Order = CblasRowMajor;
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
	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
	const enum CBLAS_ORDER Order = CblasRowMajor;
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
	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const enum CBLAS_ORDER Order = CblasRowMajor;

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
	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const enum CBLAS_ORDER Order = CblasRowMajor;

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

}
