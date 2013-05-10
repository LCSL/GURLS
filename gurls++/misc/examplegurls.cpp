/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-1013, IIT@MIT Lab
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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <ctime>
#include <cassert>
#include <fstream>

#include "gurls.h"
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


using namespace gurls;
using namespace std;

//typedef float T;
typedef double T;


int main(int argc, char *argv[])
{
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


    try
    {
        gMat2D<T> Xtr, Xte, ytr, yte;

        Xtr.readCSV(input_folder+"/Xtr.txt");
        Xte.readCSV(input_folder+"/Xte.txt");
        ytr.readCSV(input_folder+"/ytr.txt");
        yte.readCSV(input_folder+"/yte.txt");

//        cout << "Xtr = " << endl << Xtr << endl;
//        cout << "Xte = " << endl << Xte << endl;


        GurlsOptionsList opt("ExampleExperiment", true);

        OptTaskSequence *seq = new OptTaskSequence();
        *seq << "paramsel:fixsiglam" << "kernel:rbf" << "optimizer:rlsdual" << "predkernel:traintest" << "pred:dual" << "perf:macroavg";
        opt.addOpt("seq", seq);


        GurlsOptionsList * process = new GurlsOptionsList("processes", false);

        OptProcess* process1 = new OptProcess();
        *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave << GURLS::ignore << GURLS::ignore << GURLS::ignore;
        process->addOpt("one", process1);

        OptProcess* process2 = new OptProcess();
        *process2 << GURLS::load << GURLS::load << GURLS::load << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
        process->addOpt("two", process2);

        opt.addOpt("processes", process);

        std::cout << opt << std::endl;

        std::string jobid1("one");
        std::string jobid2("two");

        GURLS G;

        G.run(Xtr, ytr, opt, jobid1);

        std::cout << std::endl;

        G.run(Xte, yte, opt, jobid2);

        std::cout << opt << std::endl;

//        opt.save("par1.txt");
//        GurlsOptionsList *s1 = new GurlsOptionsList("dummy", false);
//        s1->load("par1.txt");
//        std::cout << *s1 << std::endl;

    }
    catch (gException& e)
    {
        cout << e.getMessage() << endl;
        return -1;
    }


    return 0;

}

