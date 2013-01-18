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


#ifndef _GURLS_BIGMACROAVG_H_
#define _GURLS_BIGMACROAVG_H_

#include "optlist.h"
#include "bigperf.h"
#include "utils.h"
#include "optmatrix.h"
#include "bigarray.h"
#include "bigmath.h"

namespace gurls
{

/**
 * \ingroup Performance
 * \brief BigPerfMacroAvg is the sub-class of BigPerformance that evaluates prediction accuracy
 */

template <typename T>
class BigPerfMacroAvg: public BigPerformance<T>
{

public:
    /**
     * Evaluates the average accuracy per class.
     *
     * \param X input data bigarray
     * \param Y labels bigarray
     * \param opt options with the following:
     *  - pred (settable with the class Prediction and its subclasses)
     *  - nb_pred
     *
     * This task parallelizes the costs computation
     *
     * \return a GurslOptionList equal to the field pred of opt, with the following fields added or substituted:
     *  - acc = array of prediction accuracy for each class
     *  - forho = acc
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList* BigPerfMacroAvg<T>::execute(const BigArray<T>& /*X*/, const BigArray<T>& Y, const GurlsOptionsList& opt) throw(gException)
{
    GurlsOptionsList* perf;

    if(opt.hasOpt("perf"))
    {
        GurlsOptionsList* tmp = new GurlsOptionsList("tmp");

        tmp->copyOpt("perf", opt);

        perf = tmp->getOptAs<GurlsOptionsList>("perf");

        tmp->removeOpt("perf", false);
        delete tmp;

        perf->removeOpt("acc");
        perf->removeOpt("forho");
    }
    else
        perf = new GurlsOptionsList("perf");

//    nb_pred = opt.nb_pred;
    const unsigned long nb_pred = static_cast<unsigned long>(opt.getOptAsNumber("nb_pred"));

    const BigArray<T>& pred = opt.getOptValue<OptMatrix<BigArray<T> > >("pred.pred");


//    T = y.Sizes();
//    T = T{1};
    const unsigned long t = Y.cols();

//    n_class = zeros(1,T);
    T* nClass = new T[t];
    set(nClass, (T)0.0, t);

//    flatcost = zeros(1,T);
    T* flatCost = new T[t];
    set(flatCost, (T)0.0, t);


    int numprocs;
    int myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    unsigned long blockSize = Y.rows()/numprocs;
    unsigned long remainder = 0;
    if(myid == numprocs-1)
        remainder = Y.rows()%numprocs;

    const unsigned long block_rows = blockSize + remainder;


    gMat2D<T> *y_block = new gMat2D<T>(block_rows, t);
    gMat2D<T> *ypred_block = new gMat2D<T>(block_rows, t);

//        y_block = y(i1 : i2, : );
    Y.getMatrix(myid*blockSize, 0, *y_block);

//        ypred_block = opt.pred(i1 : i2, : );
    pred.getMatrix(myid*blockSize, 0, *ypred_block);

    T* work = new T[y_block->getSize()];

//        [dummy,IY]= max(y_block,[],2);
    unsigned long* IY = new unsigned long[y_block->rows()];
    indicesOfMax(y_block->getData(), y_block->rows(), y_block->cols(), IY, work, 2);
    delete[] work;

//        [dummy,IYpred]= sort(ypred_block,2,'descend');
    unsigned long* IYpred = new unsigned long[y_block->getSize()];
    delete y_block;


    T* values = NULL;
    sort(ypred_block->getData(), ypred_block->rows(), ypred_block->cols(), &gt, values, IYpred);
    delete ypred_block;

//        for i=1:size(y_block,1)
    for(unsigned long i=0; i< block_rows; ++i)
    {
        const unsigned long index = IY[i];

//            flatcost(IY(i)) = flatcost(IY(i)) + ~ismember(IY(i),IYpred(i,1:nb_pred));
        const unsigned long value = IY[i];
        for(unsigned long j=0; j<nb_pred; ++j)
        {
            if(IYpred[i+(t*j)] == value)
            {
                ++flatCost[index];
                break;
            }
        }

//            n_class(IY(i)) = n_class(IY(i)) + 1;
        ++nClass[index];
    }

    delete[] IYpred;


    T* allNClass = new T[t];
    MPI_AllReduceT(nClass, allNClass, t, MPI_SUM, MPI_COMM_WORLD);
    delete[] nClass;

    T* allFlatCost = new T[t];
    MPI_AllReduceT(flatCost, allFlatCost, t, MPI_SUM, MPI_COMM_WORLD);
    delete[] flatCost;

//    p.acc = 1-(flatcost./n_class);
    gMat2D<T>* acc_mat = new gMat2D<T>(1, t);
    set(acc_mat->getData(), (T)1.0, t);

    rdivide(allFlatCost, allNClass, allFlatCost, t);
    delete[] allNClass;


    axpy(t, (T)-1.0, allFlatCost, 1, acc_mat->getData(), 1);
    delete[] allFlatCost;


    OptMatrix<gMat2D<T> >* acc_opt = new OptMatrix<gMat2D<T> >(*acc_mat);
    perf->addOpt("acc", acc_opt);


//    p.forho = 1-(flatcost./n_class);
    OptMatrix<gMat2D<T> >* forho_opt = new OptMatrix<gMat2D<T> >(*(new gMat2D<T>(*acc_mat)));
    perf->addOpt("forho", forho_opt);

    return perf;
}

}

#endif //_GURLS_BIGMACROAVG_H_
