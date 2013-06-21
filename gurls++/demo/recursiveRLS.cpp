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

#include "wrapper.h"
#include "rlsprimal.h"
#include "primal.h"

using namespace gurls;
typedef double T;

// size of first batch to be used for initialization
const unsigned long n0 = 100;

const T errCoeff = 1.0e-10;

const gMat2D<T>* standardRLS(const gMat2D<T> &Xtr_tot, const gMat2D<T> &ytr_tot, const gMat2D<T> &Xte_tot, const GurlsOptionsList &opt, T lambdasScaleFactor);
bool checkMatrices(const gMat2D<T>& m1, const gMat2D<T>& m2);

/**
  * Main function
  *
  * The data is already split into training and test set, and each set is
  * in the form of an input data matrix and a output labels vector.
  * Parameter selection and initial RLS estimation is carried out on a first subset of the training set.
  * Recursive RLS is run on the remainder of the training set, simulating online learning.
  * Finally the gurls++ testing process is run on the test set.
  */
int main(int argc, char* argv[])
{
    srand(time(NULL));

//    std::cout.precision(16);
//    std::cout.setf( std::ios::fixed, std::ios::floatfield);
//    std::cout.setf (std::cout.showpos);

    if(argc < 2 || argc > 3)
    {
        std::cout << "Usage: " << argv[0] << " <gurls++ data directory> [--retrain]" << std::endl;
        return EXIT_SUCCESS;
    }

    bool retrain = false;

    if(argc == 3)
        retrain = (std::string(argv[2]) == "--retrain");

    gMat2D<T> Xtr_tot, Xte_tot, ytr_tot;

    std::string XtrFileName = std::string(argv[1]) + "/Xtr.txt";
    std::string XteFileName = std::string(argv[1]) + "/Xte.txt";
    std::string ytrFileName = std::string(argv[1]) + "/ytr.txt";


    RecursiveRLSWrapper<T>* wrapper;
    if(retrain)
        wrapper = new RecursiveRLSRetrainWrapper<T>("recursiveRLSretrain");
    else
        wrapper = new RecursiveRLSWrapper<T>("recursiveRLS");

    try
    {
        std::cout << "Running recursive RLS with" << (retrain? std::string(""): std::string("out")) << " retraining..." << std::endl;

        // Load data files
        Xtr_tot.readCSV(XtrFileName);
        Xte_tot.readCSV(XteFileName);
        ytr_tot.readCSV(ytrFileName);


        // Load first batch for parameter selection and initialization
        const unsigned long n = Xtr_tot.rows();
        const unsigned long d = Xtr_tot.cols();
        const unsigned long t = ytr_tot.cols();
        const unsigned long nte = Xte_tot.rows();

        gMat2D<T> Xtr(n0, d);
        gMat2D<T> ytr(n0, t);

        unsigned long *rows = new unsigned long [n0];
        for(unsigned long i=0; i<n0; ++i)
            rows[i] = i;

        subMatrixFromRows(Xtr_tot.getData(), n, d, rows, n0, Xtr.getData());
        subMatrixFromRows(ytr_tot.getData(), n, t, rows, n0, ytr.getData());


        // Initialize model
        std::cout << "Training model..." << std::endl;
        wrapper->train(Xtr, ytr);


        // Update RLS estimator recursively
        std::cout << "\nUpdating estimator..." << std::endl;

        gVec<T> Xnew(d);
        gVec<T> ynew(t);
        for(unsigned long i=n0; i<n; ++i)
        {
            // Read a row from the file where the global training set is stored and update estimator

            getRow(Xtr_tot.getData(), n, d, i, Xnew.getData());
            getRow(ytr_tot.getData(), n, t, i, ynew.getData());

            // Update estimator with a new input pair
            wrapper->update(Xnew, ynew);
        }
        std::cout << "Update: " << tot.total_milliseconds() <<  std::endl;

        // Reinitialize model
        if(retrain)
        {
            std::cout << "Retraining model..." << std::endl;
            dynamic_cast<RecursiveRLSRetrainWrapper<T>*>(wrapper)->retrain();
        }

        // Test on independent test set
        std::cout << "Testing..." << std::endl;

        // Option 1: Test on entire data set
        gMat2D<T>* rec_result = wrapper->eval(Xte_tot);


        // Option 2: Test iteratively row-by-row
        for(unsigned long i=0; i<nte; ++i)
        {
            getRow(Xte_tot.getData(), nte, d, i, Xnew.getData());

            unsigned long index;
            double ret = wrapper->eval(Xnew, &index);

            // Check that solutions coincide
            if(!le(std::abs((*rec_result)(i, index) - ret), errCoeff*std::abs(ret)))
            {
                std::cout << "Iteration " << i << " : error detected" << std::endl;
                std::cout << " " << (*rec_result)(i, index) << " " << ret << std::endl;
            }
        }


        // Check against standard RLS
        std::cout << "Checking against standard RLS..." << std::endl;
        T alpha = retrain? (T)1.0 : ((T)n0/(T)n);
        const gMat2D<T>* std_result = standardRLS(Xtr_tot, ytr_tot, Xte_tot, wrapper->getOpt(), alpha);
        bool check = checkMatrices(*rec_result, *std_result);

        delete rec_result;
        delete std_result;

        if(check)
            return EXIT_SUCCESS;

        return EXIT_FAILURE;
    }
    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

}

const gMat2D<T>* standardRLS(const gMat2D<T> &Xtr_tot, const gMat2D<T> &ytr_tot, const gMat2D<T> &Xte_tot, const GurlsOptionsList &opt, T lambdasScaleFactor)
{
    GurlsOptionsList* stdOpt = new GurlsOptionsList(opt);
    stdOpt->removeOpt("optimizer");
    stdOpt->removeOpt("pred");

    GurlsOptionsList* paramsel = stdOpt->getOptAs<GurlsOptionsList>("paramsel");

    gMat2D<T> &lambdas = paramsel->getOptValue< OptMatrix<gMat2D<T> > >("lambdas");
    scal(lambdas.getSize(), lambdasScaleFactor, lambdas.getData(), 1);

    RLSPrimal<T> optimizer;
    stdOpt->addOpt("optimizer", optimizer.execute(Xtr_tot, ytr_tot, *stdOpt));

    PredPrimal<T> predprimal;
    gMat2D<T> empty_y;
    OptMatrix<gMat2D<T> >* pred = predprimal.execute(Xte_tot, empty_y, *stdOpt);

    const gMat2D<T>& result2 = pred->getValue();

    pred->detachValue();

    delete stdOpt;

    return &result2;
}

bool checkMatrices(const gMat2D<T>& m1, const gMat2D<T>& m2)
{
    if(m1.rows()!= m2.rows() || m1.cols()!= m2.cols())
    {
        std::ostringstream str;

        str << "Dimensions error: (" <<
                     m1.rows() << "x" << m1.cols() << ") vs (" <<
                     m2.rows() << "x" << m2.cols() << ");" << std::endl;

        throw gException(str.str());
    }

    const T* r1_it = m1.getData();
    const T* r2_it = m2.getData();
    const unsigned long size = m1.getSize();

    T min_diff = DBL_MAX;
    T max_diff = -1;
    int errorCount = 0;

    for(const T* r1_end = r1_it+size; r1_it < r1_end; ++r1_it, ++r2_it)
    {
        if(!le(std::abs(*r1_it - *r2_it), errCoeff*std::abs(*r2_it)))
        {
            T diff = fabs(*r1_it - *r2_it);
            //std::cerr << "Error " << *r1_it << " vs " << *r2_it << ", diff: " << fabs(*r1_it - *r2_it) << std::endl;

            if(lt (diff, min_diff))
                min_diff = diff;
            if(gt(diff, max_diff))
                max_diff = diff;

            ++errorCount;
        }
    }

    std::cout << "Errors found: " << errorCount << std::endl;
    if(errorCount > 0)
    {
        std::cout << " Min error: " << min_diff << std::endl << " Max error: " << max_diff << std::endl;
        return false;
    }

    return true;

}
