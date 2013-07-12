/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-2013, IIT@MIT Lab
 * All rights reserved.
 *
 * authors:  M. Santoro
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


#ifndef _GURLS_BIGSPLITHO_H_
#define _GURLS_BIGSPLITHO_H_


#include "bgurls++/bigsplit.h"
#include "bgurls++/bigarray.h"
#include "bgurls++/bigmath.h"
#include "gurls++/optmatrix.h"
#include "bgurls++/mpi_utils.h"

#include <boost/filesystem.hpp>

namespace gurls
{

/**
 * \ingroup Split
 * \brief BigSplitHo is the sub-class of BigSplit that splits data into one or more pairs of training and test samples.
 */

template <typename T>
class BigSplitHo: public BigSplit<T>
{
public:
    /**
     * Splits data into one or more pairs of training and test samples, to be used for cross-validation. The fraction of samples for the validation set is specified in the field hoproportion of opt, and the number of pairs is specified in the field nholdouts of opt
     * \param X
     * \param Y labels BigArray
     * \param opt options with the following field
     *   - hoproportion (default)
     *   - files list containing file names for BigArrays
     *   - memlimit maximum amount memory to be used performing matrix multiplications
     *
     * \return a list containing the following fields:
     *  - Xva
     *  - Yva
     *  - XvatXva
     *  - XvatYva
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList* BigSplitHo<T>::execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList &opt) throw(gException)
{
    T hoProportion = static_cast<T>(opt.getOptAsNumber("hoproportion"));

    int numprocs;
    int myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    const unsigned long n = X.rows();
    const unsigned long d = X.cols();
    const unsigned long t = Y.cols();


//    // distributed over examples (n)
//    unsigned long numRows = n/numprocs;
//    unsigned long remainder = (myid == numprocs-1)? (n%numprocs) : 0;

//    unsigned long start = myid*numRows;
//    unsigned long end = start + numRows + remainder;


//    T* work = new T[(numRows+remainder)*t];
//    unsigned long* inds = new unsigned long[numRows+remainder];
//    memset(inds, 0, (numRows+remainder)*sizeof(unsigned long));
//    gMat2D<T> mat(numRows+remainder, t);

//    Y.getMatrix(start, 0, mat);
//    indicesOfMax(mat.getData(), mat.rows(), mat.cols(), inds, work, 2);

//    delete [] work;

//    BigArray<unsigned long> indices(opt.getOptAsString("tmpfile"), n, 1);
//    indices.setMatrix(start, 0, inds, numRows+remainder, 1);
//    delete[] inds;

//    MPI_Barrier(MPI_COMM_WORLD);

//    // distributed over classes (t)
//    unsigned long numCols = t/numprocs;
//    remainder = (myid == numprocs-1)? (t%numprocs) : 0;

//    start = myid*numCols;
//    const unsigned long size = numCols + remainder;
//    end = start + size;

//    unsigned long* counts = new unsigned long[size];
//    memset(counts, 0, (size)*sizeof(unsigned long));

//    for(unsigned long i=0; i<n; ++i)
//    {
//        unsigned long ind = indices.getValue(i, 0);
//        if(ind >= start && ind < end)
//            ++counts[ind-start];
//    }


//    unsigned long **classes = new unsigned long*[size];
//    for(unsigned long i=0; i<size; ++i)
//        classes[i] = new unsigned long [counts[i]];

//    unsigned long* used = new unsigned long[size];
//    memset(used, 0, (size)*sizeof(unsigned long));

//    for(unsigned long i=0; i<n; ++i)
//    {
//        unsigned long ind = indices.getValue(i, 0);
//        if(ind >= start && ind < end)
//        {
//            unsigned long index = ind-start;
//            classes[index][used[index]++] = i;
//        }
//    }
//    delete[] used;


//    unsigned long* local_hoCounts = new unsigned long[t];
//    memset(local_hoCounts, 0, t*sizeof(unsigned long));

//    for(unsigned long i=0; i<size; ++i)
//        local_hoCounts[i+start] = static_cast<unsigned long>(std::ceil(counts[i]*hoProportion));

//    unsigned long* hoCounts = new unsigned long[t];
//    MPI_Allreduce(local_hoCounts, hoCounts, (int)t, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
//    delete[] local_hoCounts;


//    unsigned long rows = 0;
//    for(unsigned long* it = hoCounts, *end = hoCounts+t; it<end; ++it)
//        rows += *it;

//    BigArray<T>* Xva = new BigArray<T>(opt.getOptAsString("files.Xva_filename"), rows, d);
//    BigArray<T>* Yva = new BigArray<T>(opt.getOptAsString("files.Yva_filename"), rows, t);
//    MPI_Barrier(MPI_COMM_WORLD);


//    unsigned long row_offset = 0;
//    for(unsigned long i=0; i<start; ++i)
//        row_offset += hoCounts[i];

//    for(unsigned long i=0; i<size; ++i)
//    {
//        randperm(counts[i], classes[i], false);

//        for(unsigned long j=0; j<hoCounts[i+start]; ++j)
//        {
//            Xva->setRow(row_offset, X[classes[i][j]]);
//            Yva->setRow(row_offset, Y[classes[i][j]]);

//            ++row_offset;
//        }
//    }

//    for(unsigned long **it = classes, **end= classes+size; it<end; ++it)
//        delete[] *it;

//    delete[] classes;
//    delete[] counts;
//    delete[] hoCounts;



    // distributed over examples (n)
    const unsigned long numRows = n/numprocs;
    const unsigned long remainder = (myid == numprocs-1)? (n%numprocs) : 0;
    const unsigned long lastRows = numRows+remainder;


    const unsigned long rows = static_cast<unsigned long>(n*hoProportion);
    BigArray<T>* Xva = new BigArray<T>(opt.getOptAsString("files.Xva_filename"), rows, d);
    BigArray<T>* Yva = new BigArray<T>(opt.getOptAsString("files.Yva_filename"), rows, t);

    unsigned long* indices = new unsigned long[lastRows];

    randperm(lastRows, indices, true, 0);

    const unsigned long blockRows = numRows*hoProportion;
    const unsigned long lastBlockRows = lastRows*hoProportion;

    const unsigned long begin = myid*blockRows;


    gMat2D<T> X_block(lastRows, d);
    gMat2D<T> Y_block(lastRows, t);
    gMat2D<T> Xva_block(lastBlockRows, d);
    gMat2D<T> Yva_block(lastBlockRows, t);

    X.getMatrix(numRows*myid, 0, X_block);
    Y.getMatrix(numRows*myid, 0, Y_block);

    for(unsigned long i = 0; i < lastBlockRows; ++i)
    {
        const unsigned long row = indices[i];

        copy(Xva_block.getData()+i, X_block.getData()+row, d, lastBlockRows, lastRows);
        copy(Yva_block.getData()+i, Y_block.getData()+row, t, lastBlockRows, lastRows);
    }

    delete[] indices;

    Xva->setMatrix(begin, 0, Xva_block);
    Yva->setMatrix(begin, 0, Yva_block);



    MPI_Barrier(MPI_COMM_WORLD);

    BigArray<T>* XvatXva = matMult_AtB(*Xva, *Xva, opt.getOptAsString("files.XvatXva_filename"));
    BigArray<T>* XvatYva = matMult_AtB(*Xva, *Yva, opt.getOptAsString("files.XvatYva_filename"));


    GurlsOptionsList* split = new GurlsOptionsList("split");
    split->addOpt("Xva", new OptMatrix<BigArray<T> >(*Xva));
    split->addOpt("Yva", new OptMatrix<BigArray<T> >(*Yva));
    split->addOpt("XvatXva", new OptMatrix<BigArray<T> >(*XvatXva));
    split->addOpt("XvatYva", new OptMatrix<BigArray<T> >(*XvatYva));


    return split;
}

}

#endif //_GURLS_BIGSPLITHO_H_
