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

#include <iostream>
#include <string>

#include "gurls++/gprwrapper.h"

using namespace gurls;
typedef double T;


/**
  * Main function
  *
  */
int main(int argc, char* argv[])
{
    srand(static_cast<unsigned int>(time(NULL)));

//    std::cout.precision(16);
//    std::cout.setf( std::ios::fixed, std::ios::floatfield);
//    std::cout.setf (std::cout.showpos);

    if(argc < 2 || argc > 3)
    {
        std::cout << "Usage: " << argv[0] << " <data directory>" << std::endl;
        return EXIT_SUCCESS;
    }

    gMat2D<T> Xtr, Xte, ytr;

    std::string XtrFileName = std::string(argv[1]) + "/Xtr.txt";
    std::string XteFileName = std::string(argv[1]) + "/Xte.txt";
    std::string ytrFileName = std::string(argv[1]) + "/ytr.txt";

    std::string means_fileName("pred_gp_means.txt");
    std::string vars_fileName("pred_gp_vars.txt");

    GPRWrapper<T> wrapper("gp");

    try
    {
        // Load data files
        std::cout << "Loading data files..." << std::endl;
        Xtr.readCSV(XtrFileName);
        Xte.readCSV(XteFileName);
        ytr.readCSV(ytrFileName);


        wrapper.setSigma(1.0);
        wrapper.setParam(0.1);

        // Initialize model
        std::cout << "Training model..." << std::endl;
        wrapper.train(Xtr, ytr);

        // Test on independent test set
        std::cout << "Testing..." << std::endl;
        gMat2D<T> vars;
        gMat2D<T>* means = wrapper.eval(Xte, vars);

        std::cout << "Saving means matrix in " << means_fileName << "..." << std::endl;
        means->saveCSV(means_fileName);

        std::cout << "Saving variances matrix in " << vars_fileName << "..." << std::endl;
        vars.saveCSV(vars_fileName);

        delete means;

        return EXIT_SUCCESS;
    }
    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}
