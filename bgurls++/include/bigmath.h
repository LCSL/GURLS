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
 * \ingroup LinearAlgebra
 * \file
 * \brief
 */

#ifndef _GURLS_BIGMATH_H_
#define _GURLS_BIGMATH_H_

#include "gmath.h"
#include "bigarray.h"

#include <mpi/mpi.h>

namespace gurls
{

template<typename T>
int MPI_ReduceT(T *sendbuf, T *recvbuf, int count, MPI_Op op, int root, MPI_Comm comm);

// AB
template<typename T>
BigArray<T>* matMult_AB(const BigArray<T>& A, const BigArray<T>& B, const std::string& resultFile, const unsigned long memB /*const unsigned long memMB*/)
{

    if(A.cols() != B.rows())
        throw gException(Exception_Inconsistent_Size);


//    const unsigned long bytesInMB = 1024*1024;
//    const unsigned long cells = static_cast<unsigned long>((memMB*bytesInMB)/sizeof(T));
    const unsigned long cells = static_cast<unsigned long>(memB/sizeof(T));

    const unsigned long n = A.rows();
    const unsigned long d = A.cols();
    const unsigned long t = B.cols();


    int numprocs;
    int myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    BigArray<T>* ret;

    if(myid == 0)
    {
        ret = new BigArray<T>(resultFile, n, t);
        ret->setValue(0, 0, 0);
        ret->flush();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(myid != 0)
        ret = new BigArray<T>(resultFile);

//    T maxBlockSize = std::floor(static_cast<T>(cells)/d);
//    int numBlocks = static_cast<int>( std::ceil(static_cast<T>(std::max(n,t))/maxBlockSize) );
    int numBlocks = static_cast<int>( std::ceil(static_cast<T>(d*std::max(n,t))/cells) );


    const unsigned long blockRowsA = n/numBlocks;
    const unsigned long blockRowsB = t/numBlocks;
    const unsigned long lastBlockRowsA = blockRowsA + n%numBlocks;
    const unsigned long lastBlockRowsB = blockRowsB + t%numBlocks;

    gMat2D<T>* U = new gMat2D<T>;
    gMat2D<T>* V = new gMat2D<T>;

    for(int block=myid; block<numBlocks; block+=numprocs)
    {
        unsigned long rows, cols;

        if(block == (numBlocks-1))
            rows = lastBlockRowsA;
        else
            rows = blockRowsA;

        A.getMatrix(block*blockRowsA, 0, rows, d, *U);

        for(int cblock=0; cblock<numBlocks; ++cblock)
        {
            if(cblock == numBlocks -1)
                cols = lastBlockRowsB;
            else
                cols = blockRowsB;

            gMat2D<T> result(rows, cols);

            B.getMatrix(0, cblock*blockRowsB, d, cols, *V);

            dot(U->getData(), V->getData(), result.getData(), rows, d, d, cols, rows, cols, CblasNoTrans, CblasNoTrans, CblasColMajor);

            std::cout << result << std::endl;
            ret->setMatrix(block*blockRowsA, cblock*blockRowsB, result);
        }
    }

    delete U;
    delete V;

    return ret;

}

// A'B
template<typename T>
BigArray<T>* matMult_AtB(BigArray<T>& A, BigArray<T>& B, const std::string& resultFile, const unsigned long memB /*const unsigned long memMB*/)
{

    if(A.rows() != B.rows())
        throw gException(Exception_Inconsistent_Size);


//    const unsigned long bytesInMB = 1024*1024;
//    const unsigned long cells = static_cast<unsigned long>((memMB*bytesInMB)/sizeof(T));
    const unsigned long cells = static_cast<unsigned long>(memB/sizeof(T));

    const unsigned long n = A.rows();
    const unsigned long d = A.cols();
    const unsigned long t = B.cols();


    if(cells < (2*d*t)+d+t) // maximum allocated cells are 2*d*t (sum+reduced or sum+result)
        throw gException("Not enough memory available to complete the operation");


    int numprocs;
    int myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    T maxBlockSize = std::floor(static_cast<T>(cells - (d*t))/(d+t));
    int numBlocks = static_cast<int>( std::ceil(static_cast<T>(n)/maxBlockSize) );

    if(numprocs < numBlocks)
    {
        maxBlockSize = std::floor(static_cast<T>(cells - 2*(d*t))/(d+t));
        numBlocks = static_cast<int>( std::ceil(static_cast<T>(n)/maxBlockSize) );
    }

    const unsigned long blockRows = n/numBlocks;
    const unsigned long lastBlockRows = n%numBlocks;


    T* sum = new T[d*t];    // size: d*t
    set(sum, (T)0.0, d*t);

    gMat2D<T>* U = new gMat2D<T>;   // max size: blockRows*d
    gMat2D<T>* V = new gMat2D<T>;   // max size: blockRows*t

    for(int block=myid; block<numBlocks; block+=numprocs)
    {
        unsigned long rows;

        if(block == (numBlocks-1) && lastBlockRows != 0)
        {
            A.getMatrix(block*blockRows, 0, lastBlockRows, d, *U);
            B.getMatrix(block*blockRows, 0, lastBlockRows, t, *V);

            rows = lastBlockRows;
        }
        else
        {
            A.getMatrix(block*blockRows, 0, blockRows, d, *U);
            B.getMatrix(block*blockRows, 0, blockRows, t, *V);

            rows = blockRows;
        }


        if(numprocs < numBlocks)
        {
            T* result = new T[d*t]; // size: d*t

            dot(U->getData(), V->getData(), result, rows, d, rows, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);
            axpy(d*t, (T)1.0, result, 1, sum, 1);

            delete [] result;
        }
        else
            dot(U->getData(), V->getData(), sum, rows, d, rows, t, d, t, CblasTrans, CblasNoTrans, CblasColMajor);

    }

    delete U;
    delete V;

    T* reduced = NULL;
    if(myid == 0)
        reduced = new T[d*t]; // size: d*t

    MPI_ReduceT(sum, reduced, d*t, MPI_SUM, 0, MPI_COMM_WORLD);

    delete [] sum;


    BigArray<T>* ret;

    if(myid == 0)
    {
        ret = new BigArray<T>(resultFile, d, t);
        ret->setMatrix(0, 0, reduced, d, t);
        ret->flush();
        delete[] reduced;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(myid != 0)
        ret = new BigArray<T>(resultFile);

    return ret;

}

// AB'
template<typename T>
BigArray<T>* matMult_ABt(BigArray<T>& bU, BigArray<T>& bV,  const std::string& resultFile)
{
    if(bU.cols() != bV.cols())
        throw gException(Exception_Inconsistent_Size);


    int numprocs;
    int myid;

//     Find out the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

//     Get the process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    unsigned long blockSize = bU.cols()/(unsigned long)numprocs;
    unsigned long lastBlockSize = blockSize + (bU.cols()%(unsigned long)numprocs);

    gMat2D<T>* U = new gMat2D<T>;
    gMat2D<T>* V = new gMat2D<T>;

    if(myid != numprocs-1)
    {
        bU.getMatrix(0, myid*blockSize, bU.rows(), blockSize, *U);
        bV.getMatrix(0, myid*blockSize, bV.rows(), blockSize, *V);
    }
    else
    {
        bU.getMatrix(0, myid*blockSize, bU.rows(), lastBlockSize, *U);
        bV.getMatrix(0, myid*blockSize, bV.rows(), lastBlockSize, *V);
    }


    const unsigned long rows = U->rows();
    const unsigned long cols = V->rows();

    T* result = new T[rows*cols];
    dot(U->getData(), V->getData(), result, rows, U->cols(), cols, V->cols(), rows, cols, CblasNoTrans, CblasTrans, CblasColMajor);

    delete U;
    delete V;

    T* reduced = NULL;
    if(myid == 0)
        reduced = new T[rows*cols];

    MPI_ReduceT(result, reduced, rows*cols, MPI_SUM, 0, MPI_COMM_WORLD);

    delete[] result;


    BigArray<T>* ret;

    if(myid == 0)
    {
        ret = new BigArray<T>(resultFile, bU.rows(), bV.rows());
        ret->setMatrix(0, 0, reduced, rows, cols);
        ret->flush();
        delete[] reduced;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(myid != 0)
        ret = new BigArray<T>(resultFile);


    return ret;
}

}

#endif // _GURLS_BIGMATH_H_
