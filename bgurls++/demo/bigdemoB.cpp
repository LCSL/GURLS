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

#include <bigmath.h>

using namespace gurls;
using namespace std;

//typedef float T;
typedef double T;


int main(int argc, char *argv[])
{
    if(argc < 3)
    {
        cout << "Usage: " << argv[0] << "<input dir> <shared dir> [memory]." << endl;
        cout << "\t - <input dir> is the directory where the bio_TrainTest data files reside" << endl;
        cout << "\t - <shared dir> is a directory accessible by all processes (all processes must be able to access the same path)" << endl;
        cout << "\t - [memory] is the maximum amount of memory (in bytes) used by bGurls++ to performa a matrix-matrix multiplication" << endl;

        return 0;
    }

    srand(static_cast<unsigned int>(time(NULL)));

    cout.precision(4);
    cout.width(11);
    cout << fixed << right << endl;

    string input_directory(argv[1]);
    string shared_directory(argv[2]);

    int numprocs;
    int myid;

    MPI_Init(&argc, &argv);

//     Find out the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

//     Get the process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


//    "Load" bigarrays variables for the training and test set

    BigArray<T> Xtr("Xtr.nc", 0, 0);
    Xtr.readCSV(input_directory+"/Xtr.csv");

    BigArray<T> Xte("Xte.nc", 0, 0);
    Xte.readCSV(input_directory+"/Xte.csv");

    BigArray<T> ytr("ytr.nc", 0, 0);
    ytr.readCSV(input_directory+"/ytr.csv");

    BigArray<T> yte("yte.nc", 0, 0);
    yte.readCSV(input_directory+"/yte.csv");


    BGurlsOptionsList opt("bio_demoB", shared_directory, true);
    if(argc == 4)
    {
        opt.removeOpt("memlimit");
        opt.addOpt("memlimit", new OptNumberList(atoi(argv[3])));
    }


    OptTaskSequence *seq = new OptTaskSequence();
    *seq << "bigsplit:ho" << "bigparamsel:hoprimal" << "bigoptimizer:rlsprimal" << "bigpred:primal" << "bigperf:macroavg";
    opt.addOpt("seq", seq);


    GurlsOptionsList * process = new GurlsOptionsList("processes", false);

    OptProcess* process1 = new OptProcess();
    *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::ignore << GURLS::ignore;
    process->addOpt("one", process1);

    OptProcess* process2 = new OptProcess();
    *process2 << GURLS::load << GURLS::load << GURLS::computeNsave << GURLS::computeNsave;
    process->addOpt("two", process2);

    opt.addOpt("processes", process);

    if(myid ==0)
        std::cout << opt << std::endl;

    std::string jobid1("one");
    std::string jobid2("two");

    BGURLS G;

//    Run bgurls on the training set

//    cout << "---Training..." << endl;
    G.run(Xtr, ytr, opt, jobid1);
//    cout << "---Training done..." << endl;


//    Run bgurls on the test set

//    cout << "---Testing..." << endl;
    G.run(Xte, yte, opt, jobid2);
//    cout << "---Testing done..." << endl;

    if(myid ==0)
        std::cout << opt << std::endl;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
