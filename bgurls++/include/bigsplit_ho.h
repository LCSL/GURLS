/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, IIT@MIT Lab
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


#include "bigsplit.h"
#include "bigarray.h"
#include "bigmath.h"
#include "optmatrix.h"

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
     * \param X not used
     * \param Y labels bigarray
     * \param opt options with the following field
     *   - hoproportion (default)
     *   - nholdouts (default)
     *
     * \return a list containing the following fields:
     *  - indices = nholdoutsxn matrix, each row contains the indices of training and validation samples
     *  - lasts = nholdoutsx1 array, each row contains the number of elements of training set, which will be build taking the samples corresponding to the first lasts+1 elements of indices, the remainder indices will be used for validation.
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList* BigSplitHo<T>::execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList &opt) throw(gException)
{
    double hoProportion = opt.getOptAsNumber("hoproportion");

    int numprocs;
    int myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    BigArray<T>* Xva;
    BigArray<T>* Yva;
    BigArray<T>* XvatXva;
    BigArray<T>* XvatYva;

    if(myid == 0)
    {
        Xva = new BigArray<T>(opt.getOptAsString("files.Xva_filename"), X.rows(), X.cols()*hoProportion);
        Yva = new BigArray<T>(opt.getOptAsString("files.Yva_filename"), Y.rows(), Y.cols()*hoProportion);

        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        MPI_Barrier(MPI_COMM_WORLD);

        Xva = new BigArray<T>(opt.getOptAsString("files.Xva_filename"));
        Yva = new BigArray<T>(opt.getOptAsString("files.Yva_filename"));
    }


    T* order = new T[std::max(X.cols(), Y.cols())];


    if(myid == 0)
        randperm(X.cols(), order, true, 0);

    MPI_BcastT(order, X.cols(), 0, MPI_COMM_WORLD);

    unsigned long numCols = Xva->cols()/numprocs;
    unsigned long remainder = (myid == numprocs-1)? (Xva->cols()%numprocs) : 0;

    for(unsigned long i=myid*numCols; i< (myid*numCols+numCols+remainder); ++i)
        Xva->setColumn(i, X.getColumnn(order[i]));



    if(myid == 0)
        randperm(Y.cols(), order, true, 0);

    MPI_BcastT(order, Y.cols(), 0, MPI_COMM_WORLD);

    numCols = X.cols()/numprocs;
    remainder = (myid == numprocs-1)? (Y.cols()%numprocs) : 0;

    for(unsigned long i=myid*numCols; i< (myid*numCols+numCols+remainder); ++i)
        Yva->setColumn(i, Y.getColumnn(order[i]));



    delete[] order;


    XvatXva = matMult_AtB(*Xva, *Xva, opt.getOptAsString("files.XvatXva_filename"), opt.getOptAsNumber("memlimit"));
    XvatYva = matMult_AtB(*Xva, *Yva, opt.getOptAsString("files.XvatYva_filename"), opt.getOptAsNumber("memlimit"));


    GurlsOptionsList* split = new GurlsOptionsList("split");
    split->addOpt("Xva", new OptMatrix<BigArray<T> >(*Xva));
    split->addOpt("Yva", new OptMatrix<BigArray<T> >(*Yva));
    split->addOpt("XvatXva", new OptMatrix<BigArray<T> >(*XvatXva));
    split->addOpt("XvatYva", new OptMatrix<BigArray<T> >(*XvatYva));


    return split;
}

}

#endif //_GURLS_BIGSPLITHO_H_
