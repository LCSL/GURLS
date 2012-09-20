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

#ifndef _GURLS_QUICKANDDIRTY_H_
#define _GURLS_QUICKANDDIRTY_H_

#include "gurls.h"
#include "gmat2d.h"
#include "gmath.h"
#include "optlist.h"
#include "options.h"

#include <vector>

namespace gurls
{
    /**
      * Builds a simple pipeline and executes the training process
      */
    template<typename T>
    GurlsOptionsList* gurls_train(const gMat2D<T>& X, const gMat2D<T>& y)
    {

//        T = max(y);
        const unsigned long size = y.getSize();
        unsigned long t = static_cast<unsigned long>(*std::max_element(y.getData(), y.getData()+size));

//        codes = 2*eye(T) - 1;
//        y = codes(y,:);

        gMat2D<T> Y(size, t);
        set(Y.getData(), (T)-1.0, Y.getSize());
        const T one = static_cast<T>(1.0);

        for(unsigned long i=0; i<size; ++i)
            Y(i, static_cast<unsigned long>(y.getData()[i]-1)) = one;


        GurlsOptionsList* opt = new GurlsOptionsList("quickanddirty", true);

        OptTaskSequence *seq = new OptTaskSequence();
        *seq << "kernel:linear" << "paramsel:loocvdual" << "optimizer:rlsdual" << "pred:dual" << "perf:macroavg";

        opt->addOpt("seq", seq);

        GurlsOptionsList * process = new GurlsOptionsList("processes", false);

        OptProcess* process1 = new OptProcess();
        *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave << GURLS::ignore << GURLS::ignore;
        process->addOpt("one", process1);

        OptProcess* process2 = new OptProcess();
        *process2 << GURLS::load << GURLS::load << GURLS::load << GURLS::computeNsave << GURLS::computeNsave;
        process->addOpt("two", process2);

        opt->addOpt("processes", process);


        GURLS G;

        std::string jobid1("one");
        G.run(X, Y, *opt, jobid1);

        return opt;

    }

    /**
      * Executes the testing process of an already existing pipeline
      */
    template<typename T>
    void gurls_test(const gMat2D<T>& X, const gMat2D<T>& y, GurlsOptionsList& opt)
    {
//        T = max(y);
        const unsigned long size = y.getSize();
        unsigned long t = static_cast<unsigned long>(*std::max_element(y.getData(), y.getData()+size));


//        codes = 2*eye(T) - 1;
//        y = codes(y,:);

        gMat2D<T> Y(size, t);
        set(Y.getData(), (T)-1.0, Y.getSize());
        const T one = static_cast<T>(1.0);

        for(unsigned long i=0; i<size; ++i)
            Y(i, static_cast<unsigned long>(y.getData()[i]-1)) = one;

//        opt = gurls (X, y, opt,2);
        GURLS G;

        std::string jobid2("two");
        G.run(X, Y, opt, jobid2);

//        [~,yhat] = max(opt.pred,[],2);
        const gMat2D<T>& pred_mat = OptMatrix<gMat2D<T> >::dynacast(opt.getOpt("pred"))->getValue();

        gMat2D<T>* yhat = new gMat2D<T>(1, pred_mat.rows());
        T* work = new T[pred_mat.getSize()];
        maxValues(pred_mat.getData(), pred_mat.rows(), pred_mat.cols(), yhat->getData(), work, 2);
        delete[] work;
        opt.addOpt("yhat", new OptMatrix<gMat2D<T> >(*yhat));


//        acc = opt.perf.acc;
        GurlsOptionsList* perf = GurlsOptionsList::dynacast(opt.getOpt("perf"));
        GurlsOption *acc_opt = perf->getOpt("acc");

        const gMat2D<T>& acc_mat = OptMatrix<gMat2D<T> >::dynacast(acc_opt)->getValue();

        gMat2D<T>* acc_copy = new gMat2D<T>(acc_mat);
        opt.addOpt("acc", new OptMatrix<gMat2D<T> >(*acc_copy));

    }

}
#endif //_GURLS_QUICKANDDIRTY_H_
