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
#include "bigmath.h"

#include <boost/filesystem/path.hpp>

using namespace gurls;
using namespace std;
using namespace boost::filesystem3;

//typedef float T;
typedef double T;


int main(int argc, char *argv[])
{

    if(argc < 3)
    {
        cout << "Usage: " << argv[0] << " <input dir> <shared dir>." << endl;
        cout << "\t - <input dir> is the directory where the bio_TrainTest csv data files reside" << endl;
        cout << "\t - <shared dir> is a directory accessible by all processes (all processes must be able to access the same path)" << endl;

        return 0;
    }

    srand(static_cast<unsigned int>(time(NULL)));

    cout.precision(4);
    cout.width(11);
    cout << fixed << right << endl;

    path input_directory(argv[1]);
    path shared_directory(argv[2]);

    int numprocs;
    int myid;

    MPI_Init(&argc, &argv);

//     Find out the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

//     Get the process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);



//    "Load" bigarrays variables for the training and test set

    if(myid ==0)
        cout << "Loading Xtr..." << endl;

    BigArray<T> Xtr(path(shared_directory / "Xtr.nc").native(), 0, 0);
    Xtr.readCSV(path(input_directory / "Xtr.csv").native());

    if(myid ==0)
        cout << "Loading Xte..." << endl;

    BigArray<T> Xte(path(shared_directory / "Xte.nc").native(), 0, 0);
    Xte.readCSV(path(input_directory / "Xte.csv").native());

    if(myid ==0)
        cout << "Loading ytr..." << endl;

    BigArray<T> ytr(path(shared_directory / "ytr.nc").native(), 0, 0);
    ytr.readCSV(path(input_directory / "ytr.csv").native());

    if(myid ==0)
        cout << "Loading yte..." << endl;

    BigArray<T> yte(path(shared_directory / "yte.nc").native(), 0, 0);
    yte.readCSV(path(input_directory / "yte.csv").native());


    BGurlsOptionsList opt("bio_demoB", shared_directory.native(), true);

    // remove old experiments results
    if(myid == 0)
        boost::filesystem3::remove(path(opt.getOptAsString("savefile")));

    OptTaskSequence *seq = new OptTaskSequence();
    *seq << "bigsplit:ho" << "bigparamsel:hoprimal" << "bigoptimizer:rlsprimal" << "bigpred:primal" << "bigperf:macroavg";

    opt.addOpt("seq", seq);

    GurlsOptionsList * process = new GurlsOptionsList("processes", false);

    OptProcess* process1 = new OptProcess();
    *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave << GURLS::ignore << GURLS::ignore;
    process->addOpt("one", process1);

    OptProcess* process2 = new OptProcess();
    *process2 << GURLS::load << GURLS::load << GURLS::load << GURLS::computeNsave << GURLS::computeNsave;
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

//    if(myid ==0)
//        std::cout << opt << std::endl;

//    Run bgurls on the test set
//    cout << "---Testing..." << endl;
    G.run(Xte, yte, opt, jobid2);
//    cout << "---Testing done..." << endl;

    if(myid ==0)
    {
        const gMat2D<T>& acc = opt.getOptValue<OptMatrix<gMat2D<T> > >("perf.acc");
        const int accs = acc.getSize();

        std::cout.precision(4);

        std::cout << std::endl << "Prediction accurcay is:" << std::endl;

        std::cout << "\t";
        for(int i=1; i<= accs; ++i)
            std::cout << "Class " << i << "\t";

        std::cout << std::endl << "\t";

        for(int i=0; i< accs; ++i)
            std::cout << acc.getData()[i]*100.0 << "%\t";

        std::cout << std::endl;

    }

    BigArray<T>::releaseMPIData();
    MPI_Finalize();

    return EXIT_SUCCESS;
}
