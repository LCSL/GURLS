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

/**
 * \ingroup Tutorials
 * \file
 */


//#include <cstdio>
//#include <cstdlib>
//#include <iostream>
//#include <cstring>
//#include <ctime>
//#include <cassert>
//#include <fstream>

#include <mpi/mpi.h>

#include "bgurls.h"
#include "exceptions.h"
#include "gvec.h"
#include "gmat2d.h"
#include "gmath.h"
#include "norm.h"
#include "options.h"
#include "optlist.h"


#include "paramsel.h"
#include "optimization.h"
#include "pred.h"
#include "utils.h"

#include "bigarray.h"


using namespace gurls;
using namespace std;

//typedef float T;
typedef double T;

#include <bigmath.h>

int main(int argc, char *argv[])
{
    int numprocs;
    int myid;

    MPI_Init(&argc, &argv);

//     Find out the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

//     Get the process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    BigArray<T>* U, *V;

    if(myid == 0)
    {
        U = new BigArray<T>("U.nc", 3, 4);
        U->setValue(0, 0, 0);
        U->setValue(0, 1, 3);
        U->setValue(0, 2, 6);
        U->setValue(0, 3, 9);
        U->setValue(1, 0, 1);
        U->setValue(1, 1, 4);
        U->setValue(1, 2, 7);
        U->setValue(1, 3, 10);
        U->setValue(2, 0, 2);
        U->setValue(2, 1, 5);
        U->setValue(2, 2, 8);
        U->setValue(2, 3, 11);
        U->flush();

//        V = new BigArray<T> ("V.nc", 5, 4);
//        V->setValue(0, 0, 0);
//        V->setValue(0, 1, 3);
//        V->setValue(0, 2, 6);
//        V->setValue(0, 3, 9);
//        V->setValue(1, 0, 1);
//        V->setValue(1, 1, 4);
//        V->setValue(1, 2, 7);
//        V->setValue(1, 3, 10);
//        V->setValue(2, 0, 2);
//        V->setValue(2, 1, 5);
//        V->setValue(2, 2, 8);
//        V->setValue(2, 3, 11);
//        V->setValue(3, 0, 2);
//        V->setValue(3, 1, 5);
//        V->setValue(3, 2, 8);
//        V->setValue(3, 3, 11);
//        V->setValue(4, 0, 2);
//        V->setValue(4, 1, 5);
//        V->setValue(4, 2, 8);
//        V->setValue(4, 3, 11);
//        V->flush();

//        U = new BigArray<T>("U.nc", 6, 4);
//        U->setValue(0, 0, 0);
//        U->setValue(0, 1, 3);
//        U->setValue(0, 2, 6);
//        U->setValue(0, 3, 9);
//        U->setValue(1, 0, 1);
//        U->setValue(1, 1, 4);
//        U->setValue(1, 2, 7);
//        U->setValue(1, 3, 10);
//        U->setValue(2, 0, 2);
//        U->setValue(2, 1, 5);
//        U->setValue(2, 2, 8);
//        U->setValue(2, 3, 11);
//        U->setValue(3, 0, 2);
//        U->setValue(3, 1, 5);
//        U->setValue(3, 2, 8);
//        U->setValue(3, 3, 11);
//        U->setValue(4, 0, 2);
//        U->setValue(4, 1, 5);
//        U->setValue(4, 2, 8);
//        U->setValue(4, 3, 11);
//        U->setValue(5, 0, 1);
//        U->setValue(5, 1, 4);
//        U->setValue(5, 2, 7);
//        U->setValue(5, 3, 10);
//        U->flush();

        V = new BigArray<T> ("V.nc", 6, 3);
        V->setValue(0, 0, 0);
        V->setValue(0, 1, 3);
        V->setValue(0, 2, 6);
        V->setValue(1, 0, 1);
        V->setValue(1, 1, 4);
        V->setValue(1, 2, 7);
        V->setValue(2, 0, 2);
        V->setValue(2, 1, 5);
        V->setValue(2, 2, 8);
        V->setValue(3, 0, 2);
        V->setValue(3, 1, 5);
        V->setValue(3, 2, 8);
        V->setValue(4, 0, 2);
        V->setValue(4, 1, 5);
        V->setValue(4, 2, 8);
        V->setValue(5, 0, 1);
        V->setValue(5, 1, 4);
        V->setValue(5, 2, 7);
        V->flush();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(myid != 0)
    {
        U = new BigArray<T>("U.nc");
        V = new BigArray<T> ("V.nc");
    }

//    BigArray<T>* result = ba_product_uvt(*U, *V, "bigprod.nc");
//    BigArray<T>* result = matMult_AtB(*U, *V, "bigprod.nc", 248 /*584*/);
    BigArray<T>* result = matMult_AB(*V, *U, "bigprod.nc", 248 /*584*/);

    if(myid == 0)
    {
        cout << "numprocs: " << numprocs << endl;

        cout << "bU: " << endl << (*U) << endl;
        cout << "bV: " << endl << (*V) << endl;
        cout << "result: " << endl << (*result) << endl;


        gMat2D<T> sU(U->rows(), U->cols());
        gMat2D<T> sV(V->rows(), V->cols());
//        gMat2D<T> sresult(U->rows(), V->rows());
//        gMat2D<T> sresult(U->cols(), V->cols());
        gMat2D<T> sresult(V->rows(), U->cols());

        U->getMatrix(0, 0, sU);
        V->getMatrix(0, 0, sV);

//        dot(sU.getData(), sV.getData(), sresult.getData(), sU.rows(), sU.cols(), sV.rows(), sV.cols(), sresult.rows(), sresult.cols(),
//            CblasNoTrans, CblasTrans, CblasColMajor);

//        dot(sU.getData(), sV.getData(), sresult.getData(), sU.rows(), sU.cols(), sV.rows(), sV.cols(), sresult.rows(), sresult.cols(),
//            CblasTrans, CblasNoTrans, CblasColMajor);

        dot(sV.getData(), sU.getData(), sresult.getData(), sV.rows(), sV.cols(), sU.rows(), sU.cols(), sresult.rows(), sresult.cols(),
            CblasNoTrans, CblasNoTrans, CblasColMajor);

        cout << "sU: " << endl << sU << endl;
        cout << "sV: " << endl << sV << endl;
        cout << "sresult: " << endl << sresult << endl;
    }

    MPI_Finalize();

    return 0;



//    BigArray<T>*  a = new BigArray<T>("prova.nc", 3, 4);
//    cout << "A: " << a->rows() << " " << a->cols() << endl;

//    gMat2D<T>A(3, 4);
//    for(unsigned long i=0; i< A.getSize(); ++i)
//        A.getData()[i] = i;

//    cout << "A: " << A << endl;
//    cout << "rows: " << A.rows() << endl;
//    cout << "cols: " << A.cols() << endl;

//    a->setMatrix(0,0,A);

//    cout << "a: " << (*a) << endl;

//    gMat2D<T>B(3,4);
//    a->getMatrix(0,0,B);

//    cout << endl;
//    cout << "B: " << B << endl;
//    cout << "rows: " << B.rows() << endl;
//    cout << "cols: " << B.cols() << endl;


//    gMat2D<T>C;
//    a->getMatrix(0,0,3, 4, C);

//    cout << endl;
//    cout << "C: " << C << endl;
//    cout << "rows: " << C.rows() << endl;
//    cout << "cols: " << C.cols() << endl;


//    a->setValue(0, 0, 0);
//    a->setValue(0, 1, 3);
//    a->setValue(0, 2, 6);
//    a->setValue(0, 3, 9);
//    a->setValue(1, 0, 1);
//    a->setValue(1, 1, 4);
//    a->setValue(1, 2, 7);
//    a->setValue(1, 3, 10);
//    a->setValue(2, 0, 2);
//    a->setValue(2, 1, 5);
//    a->setValue(2, 2, 8);
//    a->setValue(2, 3, 11);

//    gMat2D<T>D;
//    a->getMatrix(0,0,3, 4, D);

//    cout << endl;
//    cout << "D: " << D << endl;
//    cout << "rows: " << D.rows() << endl;
//    cout << "cols: " << D.cols() << endl;

//    cout << "Rows:" << endl;
//    for(unsigned long i=0; i< a->rows(); ++i)
//    {
//        gVec<T> row = a->getRow(i);
//        cout << i << ": " << row;
//    }

//    cout << endl;

//    cout << "Columns:" << endl;
//    for(unsigned long i=0; i< a->cols(); ++i)
//    {
//        gVec<T> col = a->getColumnn(i);
//        cout << i << ": " << col;
//    }


//    a->setColumn(0, a->getColumnn(1));

//    cout << (*a) << endl;

//    a->setRow(1, a->getRow(2));

//    cout << (*a) << endl;
//    delete a;


//    BigArray<T>* A = new BigArray<T>("prova.nc", 3, 4);
//    A->setValue(0, 0, 0);
//    A->setValue(0, 1, 3);
//    A->setValue(0, 2, 6);
//    A->setValue(0, 3, 9);
//    A->setValue(1, 0, 1);
//    A->setValue(1, 1, 4);
//    A->setValue(1, 2, 7);
//    A->setValue(1, 3, 10);
//    A->setValue(2, 0, 2);
//    A->setValue(2, 1, 5);
//    A->setValue(2, 2, 8);
//    A->setValue(2, 3, 11);

//    cout << *A << endl;

//    GurlsOptionsList opt1("test1", false);
//    opt1.addOpt("big", new OptMatrix<BigArray<T> >(*A));

//    cout << opt1 << endl;

//    GurlsOptionsList opt2 ("test2", false);
//    opt2.copyOpt("big", opt1);

//    cout << opt2 << endl;

//    opt1.save("test1.bin");

//    GurlsOptionsList opt3 ("test3", false);
//    opt3.load("test1.bin");

//    cout << opt3 << endl;

//    return 0;

/*
    srand(static_cast<unsigned int>(time(NULL)));

    cout.precision(4);
    cout.width(11);
    cout << fixed << right << endl;


    string input_folder("../data");
    if (argc < 2)
    {
        cout << "========================================================================================"<< endl
             << " WARNING: No input folder provided. " << endl
             << " Using the default option: \'" << input_folder <<  "\'" << endl
             <<    " Please make sure such folder exists or provide a valid path using the following syntax:" << endl
             << " " << argv[0]
             << " <path-to-a-valid-input-folder>" << endl
             << "========================================================================================" << endl << endl;
    }
    else
        input_folder = argv[1];

    int numprocs;
    int myid;

    MPI_Init(&argc, &argv);

//     Find out the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

//     Get the process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

//    if(myid ==0)
//    {
//        cout << myid << endl;
//    try
//    {
        gMat2D<T> aux;

        aux.readCSV(input_folder+"/Xtr.txt");
        BigArray<T> Xtr("Xtr.nc", aux);

        aux.readCSV(input_folder+"/Xte.txt");
        BigArray<T> Xte("Xte.nc", aux);

        aux.readCSV(input_folder+"/ytr.txt");
        BigArray<T> ytr("ytr.nc", aux);

        aux.readCSV(input_folder+"/yte.txt");
        BigArray<T> yte("yte.nc", aux);


        GurlsOptionsList opt("ExampleExperiment", true);
        opt.addOpt("shared_dir", "/home/bozzo/Desktop/bgurls_shared_dir/");

        OptTaskSequence *seq = new OptTaskSequence();
        *seq << "paramsel:fixsiglam" << "kernel:rbf" << "optimizer:rlsdual" << "predkernel:traintest" << "pred:dual" << "perf:macroavg";
//        *seq << "bigparamsel:hoprimal";
        opt.addOpt("seq", seq);


        GurlsOptionsList * process = new GurlsOptionsList("processes", false);

        OptProcess* process1 = new OptProcess();
        *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave << GURLS::ignore << GURLS::ignore << GURLS::ignore;
//        *process1 << GURLS::computeNsave;
        process->addOpt("one", process1);

        OptProcess* process2 = new OptProcess();
        *process2 << GURLS::load << GURLS::load << GURLS::load << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
        process->addOpt("two", process2);

        opt.addOpt("processes", process);

//        if(myid ==0)
//            std::cout << opt << std::endl;

        std::string jobid1("one");
        std::string jobid2("two");

        BGURLS G;

        G.run(Xtr, ytr, opt, jobid1, true);
//        G.run(Xtr, ytr, opt, jobid1, false);

//        if(myid ==0)
//            std::cout << std::endl;

        G.run(Xte, yte, opt, jobid2, true);

        if(myid ==0)
            std::cout << opt << std::endl;

//        opt.save("par1.txt");
//        GurlsOptionsList *s1 = new GurlsOptionsList("dummy", false);
//        s1->load("par1.txt");
//        std::cout << *s1 << std::endl;

//    }
//    catch (gException& e)
//    {
//        cout << e.getMessage() << endl;
//        return -1;
//    }
//    }

//    return 0;
    MPI_Finalize();
    return EXIT_SUCCESS;
    /**/
}
