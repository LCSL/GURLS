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


#include <iostream>

#include "gurls.h"
#include "exceptions.h"
#include "gmat2d.h"
#include "options.h"
#include "optlist.h"

#include "test_utils.h"

using namespace gurls;
using namespace std;

//typedef float T;
typedef double T;


int main(int argc, char *argv[])
{
    string xtr_file, xte_file, ytr_file, yte_file, times, performances, opt_file;

    if(argc<7)
    {
        printHelp(argc, argv);
        return 0;
    }

    xtr_file = argv[1];
    xte_file = argv[2];
    ytr_file = argv[3];
    yte_file = argv[4];
    times = argv[5];
    performances = argv[6];

    if(argc==8)
        opt_file = argv[7];

    try
    {
        gMat2D<T> *Xtr, *Xte, *ytr, *yte;

        Xtr = readFile<T>(xtr_file);
        Xte = readFile<T>(xte_file);
        ytr = readFile<T>(ytr_file);
        yte = readFile<T>(yte_file);

//        cout << endl << "Input files:"<< endl;
//        cout << " Xtr " << argv[1] <<  " rows: "<< Xtr->rows() << " cols: " << Xtr->cols() << endl;
//        cout << " Xte " << argv[2] <<  " rows: "<< Xte->rows() << " cols: " << Xte->cols() << endl;
//        cout << " ytr " << argv[3] <<  " rows: "<< ytr->rows() << " cols: " << ytr->cols() << endl;
//        cout << " yte " << argv[4] <<  " rows: "<< yte->rows() << " cols: " << yte->cols() << endl;

//        cout << "Output files:" << endl;
//        cout << "times " << argv[5] << endl;
//        cout << "performances " << argv[6] << endl << endl;


        OptTaskSequence *seq = new OptTaskSequence();

        seq->addTask("split:homulti");
//        if(Xtr->rows() <= 1000)
            seq->addTask("paramsel:siglamho");
//        else
//            seq->addTask("paramsel:fixsiglam");
        seq->addTask("kernel:rbf");
        seq->addTask("optimizer:rlsdual");
        seq->addTask("predkernel:traintest");
        seq->addTask("pred:dual");
        seq->addTask("perf:macroavg");

        const int size = 7;
        double process1[size] = {GURLS::computeNsave, GURLS::computeNsave, GURLS::computeNsave, GURLS::computeNsave, GURLS::ignore, GURLS::ignore, GURLS::ignore};
        double process2[size] = {GURLS::load, GURLS::load, GURLS::load, GURLS::load, GURLS::computeNsave, GURLS::computeNsave, GURLS::computeNsave};

        std::vector<OptNumberList*> processes;
        processes.push_back(new OptNumberList(process1, size));
        processes.push_back(new OptNumberList(process2, size));

        std::vector<std::string> jobIds;
        jobIds.push_back("one");
        jobIds.push_back("two");


        std::string name("GURLShodualrbf");
        GurlsOptionsList* opt = runTest(name, *seq, processes, jobIds, *Xtr, *ytr, *Xte, *yte);


        saveTimes(times, *opt, "split paramsel kernel rls predkernel pred perf", processes, jobIds);
        savePerformances<T>(performances, *opt);

        if(argc==8)
            opt->save(opt_file);

    }
    catch (gException& e)
    {
        cout << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}
