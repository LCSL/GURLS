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


#ifndef _GURLS_CALIBRATESGD_H_
#define _GURLS_CALIBRATESGD_H_

#include "options.h"
#include "optlist.h"
#include "gmat2d.h"
#include "gmath.h"

#include "paramsel.h"
#include "perf.h"
#include "gurls.h"

namespace gurls {

/**
 * \ingroup ParameterSelection
 * \brief ParamselCalibrateSGD is the sub-class of ParamSelection that implements parameter selection for pegasos
 */

template <typename T>
class ParamSelCalibrateSGD: public ParamSelection<T>{

public:
    /**
     * Performs parameter selection when one wants to solve the problem using rls_pegasos.
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - subsize (default)
     *  - calibfile (default)
     *  - hoperf (default)
     *  - singlelambda (default)
     *  - nlambda (default)
     *
     * \return adds the following fields to opt:
     *  - lambdas = array of values of the regularization parameter lambda minimizing the validation error for each class
     *  - W = RLS coefficient vector
     */
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void ParamSelCalibrateSGD<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{
    gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
    X_OMR.transpose(X);

    gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);


//    n_estimates = 1;
    const unsigned long n_estimates = 1;

//    [n,d] = size(X);
    const int n = X_OMR.rows();
    const int t = X_OMR.cols();

    GurlsOptionsList* tmp = new GurlsOptionsList("tmp", true);

    OptTaskSequence *seq = new OptTaskSequence();
    seq->addTask("split:ho");
    seq->addTask("kernel:linear");
    seq->addTask("paramsel:hodual");
    seq->addTask("optimizer:rlsdual");
    tmp->addOpt("seq", seq);

    tmp->addOpt("hoperf", opt.getOptAsString("hoperf"));
    tmp->addOpt("singlelambda", new OptFunction(static_cast<OptFunction*>(opt.getOpt("singlelambda"))->getName()));


    GurlsOptionsList * process = new GurlsOptionsList("processes", false);

    std::vector<double> process1;
    process1.push_back(GURLS::compute);
    process1.push_back(GURLS::compute);
    process1.push_back(GURLS::computeNsave);
    process1.push_back(GURLS::computeNsave);
    process->addOpt("one", new OptNumberList(process1));

    tmp->addOpt("processes", process);

    GURLS g;

//        sub_size = opt.subsize;
    const int subsize = static_cast<int>(opt.getOptAsNumber("subsize"));

    unsigned long* idx = new unsigned long[n]; //will use only the first subsize elements
    T* lambdas = new T[n_estimates];
    T* work = new T[std::max(n_estimates+1, std::max(X_OMR.cols(), Y_OMR.cols()))];

    gMat2D<T> Mx_tmp(t, subsize);
    gMat2D<T> My_tmp(Y_OMR.cols(), subsize);

    gMat2D<T> Mx(subsize, t);
    gMat2D<T> My(subsize, Y_OMR.cols());

    //    for i = 1:n_estimates,
    for(unsigned long i=0; i<n_estimates; ++i)
    {
//        idx = randsample(n, sub_size);
        randperm(n, idx);

//        M = X(idx,:);
        subMatrixFromRows(X.getData(), n, t, idx, subsize, Mx_tmp.getData()/*, work*/);

//        if ~exist([opt.calibfile '.mat'],'file')
//            fprintf('\n\tCalibrating...');
//            %% Step 1 : Hold out parameter selection in the dual
//            name = opt.calibfile;
//            tmp.hoperf = opt.hoperf;
//            tmp = defopt(name);
//            tmp.seq = {'split:ho','kernel:linear','paramsel:hodual','rls:dual'};
//            tmp.process{1} = [2,2,2,2];
//            tmp.singlelambda = opt.singlelambda;

//            gurls(M,y(idx,:),tmp,1);
        subMatrixFromRows(Y.getData(), Y_OMR.rows(), Y_OMR.cols(), idx, subsize, My_tmp.getData()/*, work*/);

        Mx_tmp.transpose(Mx);
        My_tmp.transpose(My);

        g.run(Mx, My, *tmp, "one");

//        end
//        fprintf('\n\tLoading existing calibration');
//        load([opt.calibfile '.mat']);
//        lambdas(i) = opt.singlelambda(opt.paramsel.lambdas);

        GurlsOptionsList* paramsel = GurlsOptionsList::dynacast(tmp->getOpt("paramsel"));
        std::vector<double>& ll = OptNumberList::dynacast(paramsel->getOpt("lambdas"))->getValue();
        T lambda = static_cast<T>((OptFunction::dynacast(opt.getOpt("singlelambda")))->getValue(&(*(ll.begin())), ll.size()));
        lambdas[i] = lambda;

//        % Add rescaling
//    end
    }

    delete[] idx;

    GurlsOptionsList* paramsel;// = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));
    if(!opt.hasOpt("paramsel"))
    {
        paramsel = new GurlsOptionsList("paramsel");
        opt.addOpt("paramsel",paramsel);
    }
    else
    {
        paramsel = GurlsOptionsList::dynacast(opt.getOpt("paramsel"));
        paramsel->removeOpt("lambdas");
        paramsel->removeOpt("W");
    }


    GurlsOptionsList* rls = GurlsOptionsList::dynacast(tmp->getOpt("optimizer"));
//    params.lambdas = mean(lambdas);
    OptNumberList* LAMBDA = new OptNumberList();
    LAMBDA->add(static_cast<double>(sumv(lambdas, n_estimates, work)/n_estimates));
    paramsel->addOpt("lambdas", LAMBDA);

//    params.W = opt.rls.W;

//    OptMatrix<gMat2D<T> >* W = static_cast<OptMatrix<gMat2D<T> >*>(rls->getOpt("W"));
//    paramsel->addOpt("W", new OptMatrix<gMat2D<T> >(W->getValue()));

    GurlsOption* W = rls->getOpt("W");
    rls->removeOpt("W", false);
    paramsel->addOpt("W", W);

    delete tmp;
    delete[] lambdas;
    delete[] work;
}


}

#endif // _GURLS_CALIBRATESGD_H_
