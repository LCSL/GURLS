#ifndef BLAS_LAPACK_H
#define BLAS_LAPACK_H


/**
 * \ingroup LinearAlgebra
 * \file
 * \brief Contains prototypes for BLAS level 1, 2, 3 and Lapack routines.
 */

namespace gurls
{
    /**
      * \enum CBLAS_ORDER
      * Matrix Order (row major or column major)
      */
    enum CBLAS_ORDER     {CblasRowMajor=9, CblasColMajor=10};

    /**
      * \enum CBLAS_TRANSPOSE
      * Transposition options (no transpose, transpose, conjugate transpose)
      */
    enum CBLAS_TRANSPOSE {CblasNoTrans=0, CblasTrans=1, CblasConjTrans=2};

    /**
      * \enum CBLAS_UPLO
      * Upper/lower options (upper, lower)
      */
    enum CBLAS_UPLO      {CblasUpper=3, CblasLower=4};

    /**
      * \enum CBLAS_DIAG
      * Diagonal options (unit, non unit)
      */
    enum CBLAS_DIAG      {CblasNonUnit=5, CblasUnit=6};

    /**
      * \enum CBLAS_SIDE
      * Side options (left, right)
      */
    enum CBLAS_SIDE      {CblasLeft=7, CblasRight=8};

    /**
     * \brief BlasUtils is a convenience class to interface with Blas
     */
    class BlasUtils
    {
    public:

        /**
         * Converts a Blas enumeration into the corresponding character used as parameter in Blas fortran routines
         */
        static char charValue(int value)
        {
            static char chars[] = "NTCULNULR";

            return chars[value];
        }

    };

}

extern "C"
{


// ------ BLAS


/**
  * \brief Prototype for Blas SDOT
  *
  * Dot product of two single precision vectors
  */
float sdot_(int *n, float *x, int *incx, float *y, int *incy);

/**
  * \brief Prototype for Blas SGEMM
  *
  * Performs one of the matrix-matrix operations:
  * \f[ C = \alpha op(A) op(B) + \beta C,\f]
  * where  \f$op(X)\f$ is one of
  * \f$op(X) = X\f$ or \f$op(X) = X^T\f$,
  * \f$\alpha\f$ and \f$\beta\f$ are scalars, and \f$A\f$, \f$B\f$ and \f$C\f$ are matrices,
  * with \f$op(A)\f$ an m by k matrix, \f$op(B)\f$ a k by n matrix
  * and \f$C\f$ an m by n matrix
  */
void sgemm_(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

/**
  * \brief Prototype for Blas SGEMV
  *
  * Performs one of the matrix-vector operations
  * \f[y = \alpha A x + \beta y\f] or \f[y = \alpha A^T x + \beta y\f]
  * where \f$\alpha\f$ and \f$\beta\f$ are scalars, \f$x\f$ and \f$y\f$ are vectors and \f$A\f$ is an m by n matrix.
  */
void sgemv_(char *trans, int *m, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);

/**
  * \brief Prototype for Blas SCOPY
  *
  * Copies a vector x to a vector y
  */
void scopy_(int *n, float *x, int *incx, float *y, int *incy);

/**
  * \brief Prototype for Blas SAXPY
  *
  * Constant times a vector plus a vector
  */
void saxpy_(int *n, float *alpha, float *x, int *incx, float *y, int *incy);

/**
  * \brief Prototype for Blas SNRM2
  *
  * Returns the euclidean norm of a vector via the function name, so that \f$ SNRM2 = \sqrt{(x^T x)} \f$
  */
float snrm2_(int *n, float *x, int *incx);

/**
  * \brief Prototype for Blas SSCAL
  *
  * Scales a vector by a constant
  */
void sscal_(int *n, float *a, float *x, int *incx);

/**
  * \brief Prototype for Blas STRSM
  *
  * Solves one of the matrix equations
  * \f[op(A) X = \alpha B\f] or \f[X op(A) = \alpha B\f],
  * where \f$\alpha\f$ is a scalar, \f$X\f$ and \f$B\f$ are m by n matrices, \f$A\f$ is a unit,
  * or non-unit, upper or lower triangular matrix and \f$op(A)\f$ is one  of \f$ op(A) = A\f$
  * or \f$op(A) = A^T\f$.
  * The matrix \f$X\f$ is overwritten on \f$B\f$.
  */
void strsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb);

/**
  * \brief Prototype for Blas SSWAP
  *
  * Interchanges two vectors.
  * Uses unrolled loops for increments equal to 1.
  */
void sswap_(int *n, float *sx, int *incx, float *sy, int *incy);

/**
  * \brief Prototype for Blas DDOT
  *
  * Dot product of two double precision vectors
  */
double ddot_(int *n, double *x, int *incx, double *y, int *incy);

/**
  * \brief Prototype for Blas DGEMM
  *
  * Performs one of the matrix-matrix operations:
  * \f[ C = \alpha op(A) op(B) + \beta C,\f]
  * where  \f$op(X)\f$ is one of
  * \f$op(X) = X\f$ or \f$op(X) = X^T\f$,
  * \f$\alpha\f$ and \f$\beta\f$ are scalars, and \f$A\f$, \f$B\f$ and \f$C\f$ are matrices,
  * with \f$op(A)\f$ an m by k matrix, \f$op(B)\f$ a k by n matrix
  * and \f$C\f$ an m by n matrix
  */
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

/**
  * \brief Prototype for Blas DGEMV
  *
  * Performs one of the matrix-vector operations
  * \f[y = \alpha A x + \beta y\f] or \f[y = \alpha A^T x + \beta y\f]
  * where \f$\alpha\f$ and \f$\beta\f$ are scalars, \f$x\f$ and \f$y\f$ are vectors and \f$A\f$ is an m by n matrix.
  */
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

/**
  * \brief Prototype for Blas DCOPY
  *
  * Copies a vector x to a vector y
  */
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);

/**
  * \brief Prototype for Blas DAXPY
  *
  * Constant times a vector plus a vector
  */
void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);

/**
  * \brief Prototype for Blas DNRM2
  *
  * Returns the euclidean norm of a vector via the function name, so that \f$ SNRM2 = \sqrt{(x^T x)} \f$
  */
double dnrm2_(int *n, double *x, int *incx);

/**
  * \brief Prototype for Blas DSCAL
  *
  * Scales a vector by a constant
  */
void dscal_(int *n, double *a, double *x, int *incx);

/**
  * \brief Prototype for Blas DTRSM
  *
  * Solves one of the matrix equations
  * \f[op(A) X = \alpha B\f] or \f[X op(A) = \alpha B\f],
  * where \f$\alpha\f$ is a scalar, \f$X\f$ and \f$B\f$ are m by n matrices, \f$A\f$ is a unit,
  * or non-unit, upper or lower triangular matrix and \f$op(A)\f$ is one  of \f$ op(A) = A\f$
  * or \f$op(A) = A^T\f$.
  * The matrix \f$X\f$ is overwritten on \f$B\f$.
  */
void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);

/**
  * \brief Prototype for Blas DSWAP
  *
  * Interchanges two vectors.
  * Uses unrolled loops for increments equal to 1.
  */
void dswap_(int *n, double *sx, int *incx, double *sy, int *incy);



// ------ LAPACK



/**
  * \brief Prototype for Lapack SGETRF
  *
  * Computes an LU factorization of a general M-by-N matrix \f$A\f$
  * using partial pivoting with row interchanges.
  * The factorization has the form
  * \f[A = P L U\f]
  * where \f$P\f$ is a permutation matrix, \f$L\f$ is lower triangular with unit
  * diagonal elements (lower trapezoidal if m > n), and \f$U\f$ is upper
  * triangular (upper trapezoidal if m < n)
  */
int sgetrf_(int *m, int *n, float *a, int *lda, int *ipiv, int *info);

/**
  * \brief Prototype for Lapack SGESVD
  *
  * Computes the singular value decomposition (SVD) of a real
  * M-by-N matrix \f$A\f$, optionally computing the left and/or right singular
  * vectors. The SVD is written
  * \f[A = U \Sigma V^T \f]
  * where \f$\Sigma\f$ is an M-by-N matrix which is zero except for its
  * min(m,n) diagonal elements, \f$U\f$ is an M-by-M orthogonal matrix, and
  * V is an N-by-N orthogonal matrix. The diagonal elements of \f$\Sigma\f$
  * are the singular values of \f$A\f$; they are real and non-negative, and
  * are returned in descending order. The first min(m,n) columns of
  * \f$U\f$ and \f$V\f$ are the left and right singular vectors of \f$A\f$.
  * Note that the routine returns \f$V^T\f$, not \f$V\f$.
  */
int sgesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack DGESVD
  *
  * Computes the singular value decomposition (SVD) of a real
  * M-by-N matrix \f$A\f$, optionally computing the left and/or right singular
  * vectors. The SVD is written
  * \f[A = U \Sigma V^T \f]
  * where \f$\Sigma\f$ is an M-by-N matrix which is zero except for its
  * min(m,n) diagonal elements, \f$U\f$ is an M-by-M orthogonal matrix, and
  * V is an N-by-N orthogonal matrix. The diagonal elements of \f$\Sigma\f$
  * are the singular values of \f$A\f$; they are real and non-negative, and
  * are returned in descending order. The first min(m,n) columns of
  * \f$U\f$ and \f$V\f$ are the left and right singular vectors of \f$A\f$.
  * Note that the routine returns \f$V^T\f$, not \f$V\f$.
  */
int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack SGEEV
  *
  * Computes for an N-by-N real nonsymmetric matrix A, the
  * eigenvalues and, optionally, the left and/or right eigenvectors.
  * The right eigenvector v(j) of A satisfies
  * \f[A v(j) = \lambda (j) v(j)\f]
  * where \f$\lambda (j)\f$ is its eigenvalue.
  * The left eigenvector \f$u(j)\f$ of \f$A\f$ satisfies
  * \f[u(j)^T A = \lambda (j) u(j)^T\f]
  * where \f$u(j)^T\f$ denotes the transpose of \f$u(j)\f$.
  * The computed eigenvectors are normalized to have Euclidean norm
  * equal to 1 and largest component real.
  */
int sgeev_(char *jobvl, char *jobvr, int *n, float *a, int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack DGEEV
  *
  * Computes for an N-by-N real nonsymmetric matrix A, the
  * eigenvalues and, optionally, the left and/or right eigenvectors.
  * The right eigenvector v(j) of A satisfies
  * \f[A v(j) = \lambda (j) v(j)\f]
  * where \f$\lambda (j)\f$ is its eigenvalue.
  * The left eigenvector \f$u(j)\f$ of \f$A\f$ satisfies
  * \f[u(j)^T A = \lambda (j) u(j)^T\f]
  * where \f$u(j)^T\f$ denotes the transpose of \f$u(j)\f$.
  * The computed eigenvectors are normalized to have Euclidean norm
  * equal to 1 and largest component real.
  */
int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack SGETRI
  *
  * Computes the inverse of a matrix using the LU factorization
  * computed by SGETRF.
  *
  * This method inverts \f$U\f$ and then computes \f$A^{-1}\f$ by solving the system
  * \f$A^{-1} L = U^{-1}\f$ for \f$A^{-1}\f$.
  */
int sgetri_(int *n, float *a, int* lda, int *ipiv, float*work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack SPOTRF
  *
  * Computes the Cholesky factorization of a real symmetric positive definite matrix A.
  *
  * The factorization has the form
  * \f[A = U^T U\f],  if UPLO = 'U',
  * or
  * \f[A = L  L^T\f],  if UPLO = 'L',
  * where \f$U\f$ is an upper triangular matrix and \f$L\f$ is lower triangular
  */
int spotrf_(char *UPLO, int *n, float *a, int *lda , int *info);

/**
  * \brief Prototype for Lapack DPOTRF
  *
  * Computes the Cholesky factorization of a real symmetric positive definite matrix A.
  *
  * The factorization has the form
  * \f[A = U^T U\f],  if UPLO = 'U',
  * or
  * \f[A = L  L^T\f],  if UPLO = 'L',
  * where \f$U\f$ is an upper triangular matrix and \f$L\f$ is lower triangular
  */
int dpotrf_(char *UPLO, int *n, double *a, int *lda , int *info);

/**
  * \brief Prototype for Lapack SGELSS
  *
  * Computes the minimum norm solution to a real linear least squares problem:
  *
  *  Minimize \f$||(| b - Ax |)||_2\f$.
  *
  * Using the singular value decomposition (SVD) of \f$A\f$. \f$A\f$ is an M-by-N
  * matrix which may be rank-deficient.
  *
  * Several right hand side vectors b and solution vectors x can be
  * handled in a single call; they are stored as the columns of the
  * M-by-NRHS right hand side matrix \f$B\f$ and the N-by-NRHS solution matrix
  * \f$X\f$.
  *
  * The effective rank of \f$A\f$ is determined by treating as zero those
  * singular values which are less than RCOND times the largest singular
  * value.
  */
int sgelss_( int *m, int *n, int* nrhs, float *a, int *lda, float* b, int *ldb, float *s, float *rcond, int *rank, float *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack DGELSS
  *
  * Computes the minimum norm solution to a real linear least squares problem:
  *
  *  Minimize \f$||(| b - Ax |)||_2\f$.
  *
  * Using the singular value decomposition (SVD) of \f$A\f$. \f$A\f$ is an M-by-N
  * matrix which may be rank-deficient.
  *
  * Several right hand side vectors b and solution vectors x can be
  * handled in a single call; they are stored as the columns of the
  * M-by-NRHS right hand side matrix \f$B\f$ and the N-by-NRHS solution matrix
  * \f$X\f$.
  *
  * The effective rank of \f$A\f$ is determined by treating as zero those
  * singular values which are less than RCOND times the largest singular
  * value.
  */
int dgelss_( int *m, int *n, int* nrhs, double *a, int *lda, double* b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack SSYEV
  *
  * Computes all eigenvalues and, optionally, eigenvectors of a
  * real symmetric matrix A.
  */
int ssyev_( char* jobz, char* uplo, int* n, float* a, int* lda, float* w, float* work, int* lwork, int* info );

/**
  * \brief Prototype for Lapack DSYEV
  *
  * Computes all eigenvalues and, optionally, eigenvectors of a
  * real symmetric matrix A.
  */
int dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

/**
  * \brief Prototype for Lapack SGEQP3
  *
  * Computes a QR factorization with column pivoting of a matrix \f$A\f$: \f$A P = Q R\f$  using Level 3 BLAS.
  */
void sgeqp3_( int *m, int *n, float *A, int *lda, int *jpvt, float *tau, float *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack DGEQP3
  *
  * Computes a QR factorization with column pivoting of a matrix \f$A\f$: \f$A P = Q R\f$  using Level 3 BLAS.
  */
void dgeqp3_( int *m, int *n, double *A, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack SORGQR
  *
  * Generates an M-by-N real matrix \f$Q\f$ with orthonormal columns,
  * which is defined as the first N columns of a product of K elementary
  * reflectors of order M
  *
  *   \f[Q  =  H(1) H(2) . . . H(k)\f]
  *
  * as returned by SGEQRF.
  */
void sorgqr_(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);

/**
  * \brief Prototype for Lapack DORGQR
  *
  * Generates an M-by-N real matrix \f$Q\f$ with orthonormal columns,
  * which is defined as the first N columns of a product of K elementary
  * reflectors of order M
  *
  *   \f[Q  =  H(1) H(2) . . . H(k)\f]
  *
  * as returned by SGEQRF.
  */
void dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

}

#include"gurls++/blas_lapack.hpp"

#endif //BLAS_LAPACK_H
