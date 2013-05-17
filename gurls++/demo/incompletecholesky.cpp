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

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include "wrapper.h"
#include "optmatrix.h"
#include "icholwrapper.h"

using namespace gurls;
typedef double T;

// size of first batch to be used for initialization
const unsigned long n0 = 50;



int main(int argc, char* argv[])
{
    srand(time(NULL));

//    std::cout.precision(15);
//    std::cout.setf( std::ios::fixed, std::ios::floatfield);
//    std::cout.setf (std::cout.showpos);

    if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <gurls++ data directory>" << std::endl;
        return EXIT_SUCCESS;
    }

    gMat2D<T> Xtr_tot, Xte_tot, ytr_tot;

    std::string XtrFileName = std::string(argv[1]) + "/Xtr.txt";
    std::string XteFileName = std::string(argv[1]) + "/Xte.txt";
    std::string ytrFileName = std::string(argv[1]) + "/ytr.txt";

    ICholWrapper wrapper("ichol");
    wrapper.setSigma(13);
    wrapper.setRank(35);

    try
    {
        std::cout << "Running incomplete cholesky..." << std::endl;

        // Load data files
        Xtr_tot.readCSV(XtrFileName);
        Xte_tot.readCSV(XteFileName);
        ytr_tot.readCSV(ytrFileName);


        // Load first batch for parameter selection and initialization
        const unsigned long n = Xtr_tot.rows();
        const unsigned long d = Xtr_tot.cols();
        const unsigned long t = ytr_tot.cols();
//        const unsigned long nte = Xte_tot.rows();

        gMat2D<T> Xtr(n0, d);
        gMat2D<T> ytr(n0, t);

        unsigned long *rows = new unsigned long [n0];
        for(unsigned long i=0; i<n0; ++i)
            rows[i] = i;

        subMatrixFromRows(Xtr_tot.getData(), n, d, rows, n0, Xtr.getData());
        subMatrixFromRows(ytr_tot.getData(), n, t, rows, n0, ytr.getData());


        // Initialize model
        std::cout << "Training model..." << std::endl;
        wrapper.train(Xtr, ytr);

        std::cout << "Alpha: " <<std::endl;
        std::cout << wrapper.getOpt().getOptValue<OptMatrix<gMat2D<T> > >("paramsel.alpha") << std::endl;

//        // Test on independent test set
//        std::cout << "Testing..." << std::endl;

//        // Option 1: Test on entire data set
//        gMat2D<T>* rec_result = wrapper.eval(Xte_tot);


//        // Option 2: Test iteratively row-by-row
//        for(unsigned long i=0; i<nte; ++i)
//        {
//            getRow(Xte_tot.getData(), nte, d, i, Xnew.getData());

//            unsigned long index;
//            double ret = wrapper.eval(Xnew, &index);

//            // Check that solutions coincide
//            if(!le(std::abs((*rec_result)(i, index) - ret), errCoeff*std::abs(ret)))
//            {
//                std::cout << "Iteration " << i << " : error detected" << std::endl;
//                std::cout << " " << (*rec_result)(i, index) << " " << ret << std::endl;
//            }
//        }

//        delete rec_result;

        return EXIT_SUCCESS;
    }
    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

}
