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
 * \ingroup Tutorials
 * \file
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <ctime>
#include <cassert>

#include "gvec.h"
#include "gmat2d.h"
#include "gmath.h"
#include "norm.h"
#include "options.h"
#include "optlist.h"


#include "loocvprimal.h"
#include "rlsprimal.h"
#include "pred.h"
#include "primal.h"
#include "utils.h"

using namespace gurls;
using namespace std;

// Change this line to make a test with different data types
//typedef float T;
typedef float T;
const int Ncols = 12;
const int Nrows = 6;
const int Nsquared = 10;
const T tolerance = static_cast<T>(1e-4);


int main(int argc, char *argv[])
{

    cout.precision(4);
    cout.width(11);
    cout << fixed << right << endl;

    try{
        gMat2D<T> A(Nrows,Ncols);
        A = 2;
        cout << "A = " << endl << A << endl;

        gMat2D<T> B(A.cols(), A.rows());
        A.transpose(B);
        cout << "A\' = " << endl <<  B << endl;

        gVec<T> X(A.cols());
        X = 1;
        gVec<T> Y(A.rows());

        cout << "X = " << X << endl;
        cout << "Y = " << Y << endl;

        cout << "L0 norm of Y = " << norm(Y, L0norm) << endl << endl;
        cout << "L1 norm of X = " << norm(X, L1norm) << endl << endl;
        cout << "L2 norm of X = " << norm(X, L2norm) << endl << endl;

        dot(A,X,Y);
        cout << "Y = " << Y << endl;

        B = 3;
        cout << "B = " << endl << B << endl;
        dot(B,Y,X);
        cout << "B*Y = " << X << endl;

        gMat2D<T> C(A.rows(), A.rows());
        dot(A,B,C);

        cout << "A*B = " << endl << C << endl;
        cout << "A.*B " << endl << (A*B) << endl;

        T thres = static_cast<T>(0.81);

        B.randomize();
        cout << "randomize(B) = " << endl << B << endl;

        cout << "sum(B) = " << endl << B.sum() << endl;
        cout << "B<"<< thres<<") = " << endl << (B<thres) << endl;
        cout << (B<thres).sum() << endl;
        cout << "B>="<< thres<<") = " << endl << (B>=thres) << endl;
        cout << (B>=thres).sum() << endl;


        cout << "B.where(B<" << thres << ") = " << endl << B.where(B < thres).asMatrix() << endl;
        cout << "B.where(B>=" << thres << ") = " << endl << B.where(B >= thres).asMatrix() << endl;
        cout << "min(B.where(B>=" << thres << ")) = " << endl << B.where(B >= thres).min() << endl;
        cout << "max(B.where(B>=" << thres << ")) = " << endl << B.where(B >= thres).max() << endl;
        cout << "min(B.where(B<" << thres << ")) = " << endl << B.where(B < thres).min() << endl;
        cout << "max(B.where(B<" << thres << ")) = " << endl << B.where(B < thres).max() << endl;

        cout << endl << endl;

        cout << "B.where(B>=" << thres << ") = " << endl << B.where(B.compare(thres,gMat2D<T>::GreaterEq)) << endl;
        cout << "min(B.where(B>=" << thres << ") = " << endl << (B.where(B.compare(thres,gMat2D<T>::GreaterEq))).min() << endl;
        cout << "max(B.where(B>=" << thres << ") = " << endl << (B.where(B.compare(thres,gMat2D<T>::GreaterEq))).max() << endl;

        cout << endl << endl;

        gVec<T> boolOut = B.where(B >= thres);
        cout << "B >= " << thres << endl << boolOut << endl;

        cout << endl << endl;
        cout << "B = " << endl << B << endl;
        cout << "B(2,:) = " << endl << B[2].asMatrix(false) << endl;
        cout << "B(:,2) = " << endl << B(2).asMatrix() << endl;


        cout << endl << endl;
        cout << "B = " << endl << B << endl;
        cout << "TESTING B(2,:) = 2" << endl;
        X = 2;
        B.setRow(X, 2);
        cout << "B(2,:) = " << endl << B[2].asMatrix(false) << endl;
        cout << "B(:,2) = " << endl << B(2).asMatrix() << endl;

        cout << "B = " << endl << B << endl;

        std::cout << "min(B) COLUMNWISE = " << endl << B.min(gurls::COLUMNWISE)->asMatrix(false) << endl;
        std::cout << "argmin(B) COLUMNWISE = " << endl << B.argmin(gurls::COLUMNWISE)->asMatrix(false) << endl;
        std::cout << "min(B) ROWWISE = " << endl << B.min(gurls::ROWWISE)->asMatrix(true) << endl;
        std::cout << "argmin(B) ROWWISE = " << endl << B.argmin(gurls::ROWWISE)->asMatrix(true) << endl;
        std::cout << "max(B) COLUMNWISE = " << endl << B.max(gurls::COLUMNWISE)->asMatrix(false) << endl;
        std::cout << "argmax(B) COLUMNWISE = " << endl << B.argmax(gurls::COLUMNWISE)->asMatrix(false) << endl;
        std::cout << "max(B) ROWWISE = " << endl << B.max(gurls::ROWWISE)->asMatrix(true) << endl;
        std::cout << "argmax(B) ROWWISE = " << endl << B.argmax(gurls::ROWWISE)->asMatrix(true) << endl;
        std::cout << "sum(B), COLUMNWISE = " << endl << B.sum(gurls::COLUMNWISE)->asMatrix(false) << endl;
        std::cout << "sum(B) ROWWISE = " << endl << B.sum(gurls::ROWWISE)->asMatrix(true) << endl;

        assert(((B.compare(thres,gMat2D<T>::Less))+(B.compare(thres,gMat2D<T>::GreaterEq))).allEqualsTo(true));

        cout << "1./B = " << endl << B.reciprocal() << endl;

        cout << endl
             << "====================================================" << endl
             << " Computing the pseusodinverse of a general matrix M" << endl
             << "====================================================" << endl
             << endl;

        gMat2D<T>* Binv = new gMat2D<T>(B.cols(),B.rows());
        *Binv = 0;
        pinv(B, *Binv);
        cout << "B = " << endl << B << endl;
        cout << "pinv(B) = " << endl << *Binv << endl;

        gMat2D<T> *TMP1;
        if (B.rows() >= B.cols()){

            TMP1 = new gMat2D<T>(Binv->rows(), B.cols());
            dot(*Binv, B, *TMP1);

        } else {
            TMP1 = new gMat2D<T>(B.rows(), Binv->cols());
            dot(B, *Binv, *TMP1);
        }
        cout << "pinv(B)*B = "  << endl << *TMP1 << endl;
        // WARNING: check why the program halt with a segmentation fault is the following line is de-commented
        //delete[] TMP1;

        cout << endl
             << "===================================================" << endl
             << " Performing svd decomposition of a square matrix M" << endl
             << "===================================================" << endl
             << endl;

        gMat2D<T> *M = new gMat2D<T>(Nsquared, Nsquared),
                *U = new gMat2D<T>(Nsquared, Nsquared),
                *Vt = new gMat2D<T>(Nsquared, Nsquared);

        TMP1 = new gMat2D<T>(Nsquared, Nsquared);
        gMat2D<T> *TMP2 = new gMat2D<T>(Nsquared, Nsquared);


        gVec<T> *W = new gVec<T>(Nsquared);

        M->randomize();

        cout << "M = " << endl << *M;

        cout << "[U, W, Vt] = svd (M)" << endl;
        svd(*M, *U, *W, *Vt);

        *TMP1 = 0;
        TMP1->setDiag(*W);
        cout << "U = " << endl << *U;
        cout << "W = " << endl << *TMP1 << endl;
        cout << "Vt = " << endl << *Vt << endl;

        dot(*TMP1, *Vt, *TMP2);
        dot(*U, *TMP2, *TMP1);
        cout << "U*W*Vt = " << endl << *TMP1 << endl;

        assert(M->closeTo(*TMP1, tolerance));

        delete M;
        delete U;
        delete TMP1;
        delete TMP2;
        delete Vt;
        delete W;

        cout << endl
             << "=================================================" << endl
             << " Performing cholesky decomposition of a positive " << endl
             << " definite matrix M" << endl
             << "=================================================" << endl
             << endl;

        M = new gMat2D<T>(Nsquared, Nsquared);
        U = new gMat2D<T>(Nsquared, Nsquared);
        TMP1 = new gMat2D<T>(Nsquared, Nsquared);
        TMP2 = new gMat2D<T>(Nsquared, Nsquared);

        M->randomize();
        M->transpose(*TMP1);
        dot(*M, *TMP1, *U);
        cout << "U = " << endl << *U;
        *TMP1 = 0;
        *TMP2 = 0;
        *M = 0;

        cout << "[M] = cholesky(U)" << endl;
        cholesky(*U, *M);
        cout << "M = " << endl << *M;
        M->transpose(*TMP1);
        cout << "M\' = " << endl << *TMP1;
        dot(*TMP1, *M, *TMP2);
        cout << "M\'*M = " << endl << *TMP2 << endl;

        assert(U->closeTo(*TMP2, tolerance));

        delete M;
        delete U;
        delete TMP1;
        delete TMP2;

        // STARTING TO TEST SOME MORE INTERESTING FUNCTIONALITY OF GURLS
        GurlsOptionsList opt1("LOOCV", true);

        OptStringList *l = new OptStringList();
        l->add(std::string("prova"));
        l->add(std::string("inserimento"));
        l->add(std::string("stringhe"));
        l->add(std::string("in"));
        l->add(std::string("una"));
        l->add(std::string("lista"));

        OptNumberList *nl = new OptNumberList();
        for (double i = 0; i< 15; i+=1){
            nl->add(i);
        }
        opt1.addOpt("Testing StringListOption", l);
        opt1.addOpt("Testing NumberListOption", nl);

        //double (*minus)(int,int) = subtraction
        //OptFunction *fl = new OptFunction("mean", mean);
        OptFunction *fl = new OptFunction("mean");
        opt1.addOpt("Testing FunctionOption", fl);

        cout << "=========================" << endl;
        double *prova = new double[Nsquared];
        for (int i = 0 ; i < Nsquared; i++){
            prova[i]=i+1;
        }
        cout << fl->getValue(prova, Nsquared) << endl;
        //cout << (new OptFunction("median", median))->getValue(prova, Nsquared) << endl;;
        //cout << (new OptFunction("min", min))->getValue(prova, Nsquared) << endl;
        //cout << (new OptFunction("max", max))->getValue(prova, Nsquared) << endl;
        cout << (new OptFunction("median"))->getValue(prova, Nsquared) << endl;;
        cout << (new OptFunction("min"))->getValue(prova, Nsquared) << endl;
        cout << (new OptFunction("max"))->getValue(prova, Nsquared) << endl;
        M = new gMat2D<T>(Nsquared, Nsquared);
        M->randomize();
        cout << (new OptMatrix<gMat2D<T> >(*M))->getValue() << endl;
        cout << "=========================" << endl;
        cout << "=========================" << endl;
        OptMatrix<gMat2D<T> >* OMTMP = new OptMatrix< gMat2D<T> > (*M);
        OMTMP->operator <<(cout) << endl;
        cout << "=========================" << endl;
        cout << "=========================" << endl;

        //opt1.addOpt("Testing FunctionOption (mean)", new OptFunction("mean", mean));
        //opt1.addOpt("Testing FunctionOption (median)", new OptFunction("median", median));
        //opt1.addOpt("Testing FunctionOption (max)", new OptFunction("max", max));
        //opt1.addOpt("Testing FunctionOption (min)", new OptFunction("min", min));
        opt1.addOpt("Testing FunctionOption (mean)", new OptFunction("mean"));
        opt1.addOpt("Testing FunctionOption (median)", new OptFunction("median"));
        opt1.addOpt("Testing FunctionOption (max)", new OptFunction("max"));
        opt1.addOpt("Testing FunctionOption (min)", new OptFunction("min"));

        opt1.addOpt("Testing MatrixOption", new OptMatrix<gMat2D<T> >(*M));
        opt1.addOpt("Testing MatrixOption", new OptMatrix<gMat2D<T> >(*M));

        cout << opt1 << endl;


        cout << "Testing LoocvPrimal... ";
        ParamSelLoocvPrimal<T> loocv;
        loocv.execute(A, B, opt1);
        cout << "done." << endl;

        cout << "Testing RLSPrimal... " ;
        RLSPrimal<T> rls;
        rls.execute(A, B, opt1);
        cout << "done." << endl;

        cout << "Testing PredPrimal... ";
        //PredPrimal< gMat2D <T> >pred;
        PredPrimal <T> pred;
        /*gMat2D<T>& Z = */pred.execute(*M, *M, opt1);
        cout << "done." << endl;
    }
    catch (gException& e){
        cout << e.getMessage() << endl;
        return -1;
    }

    gMat2D<T> out(1, Nsquared);
    out.randomize();
    out-=0.5;
    cout  << "output = " << endl << out << endl;
    gMat2D<bool>& gt = out>0.;
    out+=gMat2D<T>::rand(1, Nsquared);
    cout  << "noisy output = " << endl << out << endl;
    cout  << "ground truth = " << endl << gt << endl;
    cout  << "noisy classification results= " << endl << (out>0)<< endl;

    T* outPtr = out.getData();
    bool* gtPtr_bool = gt.getData();
    T *gtPtr = new T[gt.getSize()];
    for (unsigned long i = 0; i < gt.getSize(); ++i){
        gtPtr[i] = static_cast<T>(2.0)*gtPtr_bool[i]-1;
    }
    double ap = precrec_driver(outPtr, gtPtr, out.getSize());

    cout << "average precision = " << ap << endl;

    delete [] gtPtr;

    return 0;

}
