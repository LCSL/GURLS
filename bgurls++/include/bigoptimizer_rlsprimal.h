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


#ifndef _GURLS_BIGRLSPRIMAL_H_
#define _GURLS_BIGRLSPRIMAL_H_


#include "optmatrix.h"
#include "optfunction.h"
#include "bigoptimization.h"
#include "bigmath.h"
#include "utils.h"

namespace gurls
{

/**
* \ingroup Optimization
* \brief BigRLSPrimal is the sub-class of BigOptimizer that implements RLS with the primal formulation
*/
template <typename T>
class BigRLSPrimal: public BigOptimizer<T>
{

public:
   /**
    * Computes a classifier for the primal formulation of RLS.
    * The regularization parameter is set to the one found in the field paramsel of opt.
    * In case of multiclass problems, the regularizers need to be combined with the function specified inthe field singlelambda of opt
    *
    * \param X input data bigarray
    * \param Y labels bigarray
    * \param opt options with the following:
    *  - singlelambda (default)
    *  - paramsel (settable with the class ParamSelection and its subclasses)
    *  - files list containing file names for BigArrays
    *  - tmpfile path of a file used to store and load temporary data
    *  - memlimit maximum amount memory to be used performing matrix multiplications
    *
    * BigArray Multiplications are performed in parallel
    *
    * \return adds to opt the field optimizer which is a list containing the following fields:
    *  - W = matrix of coefficient vectors of rls estimator for each class
    *  - C = empty matrix
    *  - X = empty matrix
    */
   GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt);
};


template <typename T>
GurlsOptionsList* BigRLSPrimal<T>::execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList &opt)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    const std::string tempFileName = opt.getOptAsString("tmpfile");

    const BigArray<T> *bK, *bXty;

    //	K = X'*X;
    if(!opt.hasOpt("paramsel.XtX"))
        bK = matMult_AtB(X, X, opt.getOptAsString("files.XtX_filename"), opt.getOptAsNumber("memlimit"));
    else
        bK = &opt.getOptValue<OptMatrix<BigArray<T> > >("paramsel.XtX");

    //	Xty = X'*y;
    if(!opt.hasOpt("paramsel.Xty"))
        bXty = matMult_AtB(X, Y, opt.getOptAsString("files.Xty_filename"), opt.getOptAsNumber("memlimit"));
    else
        bXty = &opt.getOptValue<OptMatrix<BigArray<T> > >("paramsel.Xty");

    gMat2D<T> *W;

    if(myid == 0)
    {
        gMat2D<T> K, Xty;
        bK->getMatrix(0, 0, bK->rows(), bK->cols(), K);
        bXty->getMatrix(0, 0, bXty->rows(), bXty->cols(), Xty);

        const gMat2D<T> &ll = opt.getOptValue<OptMatrix<gMat2D<T> > >("paramsel.lambdas");
        T lambda = opt.getOptAs<OptFunction>("singlelambda")->getValue(ll.getData(), ll.getSize());

        const unsigned long n = X.rows();
        const unsigned long d = X.cols();
        const unsigned long t = Y.cols();

        W = rls_primal_driver(K.getData(), Xty.getData(), n, d, t, lambda);
        W->save(tempFileName);

        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        W = new gMat2D<T>();

        MPI_Barrier(MPI_COMM_WORLD);

        W->load(tempFileName);
    }

    if(!opt.hasOpt("paramsel.XtX"))
        delete bK;

    if(!opt.hasOpt("paramsel.Xty"))
        delete bXty;

    GurlsOptionsList* optimizer = new GurlsOptionsList("optimizer");

    optimizer->addOpt("W", new OptMatrix<gMat2D<T> >(*W));

    gMat2D<T>* emptyC = new gMat2D<T>();
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*emptyC));

    //	cfr.X = [];
    gMat2D<T>* emptyX = new gMat2D<T>();
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*emptyX));


    return optimizer;
}


}
#endif // _GURLS_BIGRLSPRIMAL_H_
