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

#include <iostream>
#include <string>

#include "gvec.h"
#include "gmat2d.h"
#include "optlist.h"
#include "quickanddirty.h"

using namespace gurls;

/**
  * Builds a matrix reading elements from a text file in CSV format Parameter \a ROWM indicates
  * whether to load the matrix in row major order or in column major one.
  */
template<typename T>
gMat2D<T> *readFile(const std::string &fileName, bool ROWM = true );

/**
  * Main function
  */
int main(int argc, char* argv[])
{
    typedef double T;

    if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <gurls++ data directory>" << std::endl;
        return EXIT_SUCCESS;
    }

    gMat2D<T> *Xtr, *Xte, *ytr, *yte;

    std::string XtrFileName = std::string(argv[1]) + "Xtr.txt";
    std::string ytrFileName = std::string(argv[1]) + "ytr_onecolumn.txt";
    std::string XteFileName = std::string(argv[1]) + "Xte.txt";
    std::string yteFileName = std::string(argv[1]) + "yte_onecolumn.txt";

    try
    {
        //load the training data
        Xtr = readFile<T>(XtrFileName);
        ytr = readFile<T>(ytrFileName);

        //load the test data
        Xte = readFile<T>(XteFileName);
        yte = readFile<T>(yteFileName);


        //train the classifer
        GurlsOptionsList* opt = gurls_train(*Xtr, *ytr);

        //predict the labels for the test set and asses prediction accuracy
        gurls_test(*Xte, *yte, *opt);


        const gMat2D<T>& acc = OptMatrix<gMat2D<T> >::dynacast(opt->getOpt("acc"))->getValue();
        const int max = static_cast<int>(*std::max_element(ytr->getData(), ytr->getData()+ytr->getSize()));
        const int accs = acc.getSize();

        std::cout.precision(4);

        std::cout << std::endl << "Prediction accurcay is:" << std::endl;

        for(int i=1; i<= max; ++i)
            std::cout << "\tClass " << i << "\t";

        std::cout << std::endl;

        for(int i=0; i< accs; ++i)
            std::cout << "\t" << acc.getData()[i]*100.0 << "%\t";

        std::cout << std::endl;

        return EXIT_SUCCESS;
    }
    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}

template<typename T>
gMat2D<T> * readFile(const std::string &fileName, bool ROWM)
{
    std::vector<std::vector< T > > matrix;
    std::ifstream in(fileName.c_str());

    int rows = 0;
    int cols = 0;
    gMat2D<T> *g;

    if(!in.is_open())
        throw gurls::gException("Cannot open file " + fileName);

    std::string line;
    while (std::getline(in, line))
    {
        std::istringstream ss(line);
        std::vector< T > tf;
        std::copy(std::istream_iterator< T >(ss), std::istream_iterator< T >(), std::back_inserter(tf));

        matrix.push_back(tf);
        ++rows;
    }
    in.close();

    cols = matrix[0].size();

    g =  new gMat2D< T >(rows, cols);
    T* buffer = g->getData();

    for(int i=0; i<rows; ++i)
    {
        for(int j=0; j<cols; ++j)
        {
            if(ROWM)
                buffer[i*cols+j]= matrix[i][j];
            else
                buffer[j*rows+i]= matrix[i][j];
        }
    }

    return g;
}
