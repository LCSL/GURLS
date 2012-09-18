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
 * \ingroup Tutorials
 * \file
 */

#include <iostream>
#include "gurls.h"
#include "exceptions.h"
#include "gmat2d.h"
#include "options.h"
#include "optlist.h"

#include "gmath.h"

using namespace gurls;
using namespace std;

typedef double T; ///< Data type of the matrices elements

/**
  * Main function
  */
int main(int argc, char *argv[])
{
    string xtr_file, xte_file, ytr_file, yte_file;

    // check that all inputs are given
    if(argc<4)
    {
        std::cout << "========================================================================================"<< std::endl
        << " Wrong parameters number ("<<argc <<")." << std::endl
        << " Provide a valid path for training, test and output files using the following syntax:" << std::endl
        << " \n\n\t " << argv[0] << " xtr xte ytr yte" << std::endl
        << "========================================================================================" << std::endl << std::endl;
        return 0;
    }

    // get file names from input
    xtr_file = argv[1];
    xte_file = argv[2];
    ytr_file = argv[3];
    yte_file = argv[4];

    try
    {
        gMat2D<T> Xtr, Xte, ytr, yte;

        // load data from file
        Xtr.readCSV(xtr_file);
        Xte.readCSV(xte_file);
        ytr.readCSV(ytr_file);
        yte.readCSV(yte_file);

        // specify the task sequence
        OptTaskSequence *seq = new OptTaskSequence();
        *seq << "paramsel:loocvprimal" << "optimizer:rlsprimal" << "pred:primal" << "perf:macroavg" << "perf:precrec";


        GurlsOptionsList * process = new GurlsOptionsList("processes", false);

        // defines instructions for training process
        OptProcess* process1 = new OptProcess();
        *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::ignore << GURLS::ignore << GURLS::ignore;
        process->addOpt("one", process1);

        // defines instructions for testing process
        OptProcess* process2 = new OptProcess();
        *process2 << GURLS::load << GURLS::load << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
        process->addOpt("two", process2);

        // build an options' structure
        GurlsOptionsList* opt = new GurlsOptionsList("GURLSlooprimal", true);
        opt->addOpt("seq", seq);
        opt->addOpt("processes", process);

        GURLS G;

        string jobId0("one");
        string jobId1("two");

        // run gurls for training
        G.run(Xtr, ytr, *opt, jobId0);

        // run gurls for testing
        G.run(Xte, yte, *opt, jobId1);

    }
    catch (gException& e)
    {
        cout << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}
