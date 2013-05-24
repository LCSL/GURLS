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
#include "splitho.h"

#include "gurls.h"

using namespace gurls;
typedef double T;


GurlsOptionsList *runPipeline(const gMat2D<T>& Xte, const gMat2D<T>&yte, const gMat2D<T> &Xtr, const gMat2D<T> &alpha, double sigma);

int main(int argc, char* argv[])
{
//    std::cout.precision(15);
//    std::cout.setf( std::ios::fixed, std::ios::floatfield);
//    std::cout.setf (std::cout.showpos);

    if(argc < 6)
    {
        std::cout << "Usage: " << argv[0] << " <data directory> <sigma> <rank_max> <n_rank> <output_dir> [--split]" << std::endl;
        return EXIT_SUCCESS;
    }

    srand(time(NULL));

    boost::posix_time::ptime begin, end;
    boost::posix_time::time_duration diff;

    gMat2D<T> Xtr_tot, Xte_tot, ytr_tot, yte_tot;

    std::string dataDir(argv[1]);
    std::string outputDir(argv[5]);

    std::string XtrFileName = dataDir + "/Xtr.txt";
    std::string XteFileName = dataDir + "/Xte.txt";
    std::string ytrFileName = dataDir + "/ytr.txt";
    std::string yteFileName = dataDir + "/yte.txt";

    bool split = false;
    if(argc == 7)
        split = (std::string(argv[6]) == "--split");

    ICholWrapper wrapper("ichol");

    try
    {
        // Load data files
        Xtr_tot.readCSV(XtrFileName);
        Xte_tot.readCSV(XteFileName);
        ytr_tot.readCSV(ytrFileName);
        yte_tot.readCSV(yteFileName);

        gMat2D<T> *Xtrain;
        gMat2D<T> *ytrain;

        if(split)
        {
            SplitHo<T> splitTask;

            GurlsOptionsList* s_opt = new GurlsOptionsList("");
            s_opt->addOpt("hoproportion", new OptNumber(0.2));
            s_opt->addOpt("nholdouts", new OptNumber(1));

            GurlsOptionsList* splitOpt = splitTask.execute(Xtr_tot, ytr_tot, *s_opt);
            const gMat2D<unsigned long>& indices = splitOpt->getOptValue<OptMatrix<gMat2D<unsigned long> > >("indices");
            const gMat2D<unsigned long>& lasts = splitOpt->getOptValue<OptMatrix<gMat2D<unsigned long> > >("lasts");

            const unsigned long ntr = lasts.getData()[0];
            const unsigned long nva = indices.getSize()-ntr;

            const unsigned long *tr = indices.getData();
            const unsigned long *va = indices.getData()+ntr;


            const unsigned long n = Xtr_tot.rows();
            const unsigned long d = Xtr_tot.cols();
            const unsigned long t = ytr_tot.cols();


            Xtrain = new gMat2D<T>(ntr, d);
            subMatrixFromRows(Xtr_tot.getData(), n, d, tr, ntr, Xtrain->getData());

            ytrain = new gMat2D<T>(ntr, t);
            subMatrixFromRows(ytr_tot.getData(), n, t, tr, ntr, ytrain->getData());

            gMat2D<T> Xva(nva, d);
            subMatrixFromRows(Xtr_tot.getData(), n, d, va, nva, Xva.getData());

            gMat2D<T> yva(nva, t);
            subMatrixFromRows(ytr_tot.getData(), n, t, va, nva, yva.getData());

            delete splitOpt;

            wrapper.setXva(Xva);
            wrapper.setyva(yva);
        }
        else
        {
            Xtrain = &Xtr_tot;
            ytrain = &ytr_tot;

            wrapper.setXva(Xte_tot);
            wrapper.setyva(yte_tot);
        }

        double sigma = atof(argv[2]);
        wrapper.setSigma(sigma);
        wrapper.setRankMax(atoi(argv[3]));
        wrapper.setNRank(atoi(argv[4]));

        gMat2D<T> times(1, 5);

        // Run 1
        std::cout << "Run 1..." << std::endl;
        begin = boost::posix_time::microsec_clock::local_time();

        wrapper.train(*Xtrain, *ytrain);

        end = boost::posix_time::microsec_clock::local_time();
        diff = end-begin;

        //std::cout << "Total time: " << diff.total_milliseconds() <<  std::endl;
        times.getData()[0] = diff.total_milliseconds();

        const gMat2D<T> &times1 = wrapper.getOpt().getOptValue<OptMatrix<gMat2D<T> > >("paramsel.times");
        times1.saveCSV(outputDir + "/time.txt");


        const gMat2D<T> &perfs1 = wrapper.getOpt().getOptValue<OptMatrix<gMat2D<T> > >("paramsel.acc");
        perfs1.saveCSV(outputDir + "/perf.txt");

        std::cout << "Time: " << times1 <<  std::endl;
        std::cout << "Perf: " << perfs1 <<  std::endl;

////        T maxPerf = wrapper.getOpt().getOptAsNumber("paramsel.maxPerf");
//        unsigned long maxRank = wrapper.getOpt().getOptAsNumber("paramsel.maxRank");


//        wrapper.setRankMax(maxRank);
//        wrapper.setNRank(1);



//        // Run 2
//        std::cout << "Run 2..." << std::endl;
//        begin = boost::posix_time::microsec_clock::local_time();
//        wrapper.train(Xtr_tot, ytr_tot);
//        end = boost::posix_time::microsec_clock::local_time();

//        diff = end-begin;

//        std::cout << "Total time: " << diff.total_milliseconds() <<  std::endl;
//        times.getData()[1] = diff.total_milliseconds();

//        const gMat2D<T> &times2 = wrapper.getOpt().getOptValue<OptMatrix<gMat2D<T> > >("paramsel.times");
//        times2.saveCSV(outputDir + "/internaltimes_2.txt");
//        std::cout << "Iteration times: " << times2 <<  std::endl;


//        std::cout << "Pipeline: " <<  std::endl;
//        const gMat2D<T> &alpha = wrapper.getOpt().getOptValue<OptMatrix<gMat2D<T> > >("paramsel.alpha");
//        GurlsOptionsList *pipeline_opt = runPipeline(Xte_tot, yte_tot, Xtr_tot, alpha, sigma);

//        gMat2D<T> &ptimes = pipeline_opt->getOptValue<OptMatrix<gMat2D<T> > >("time.one");

//        std::cout << "Times: " << ptimes <<  std::endl;
//        copy(times.getData()+2, ptimes.getData(), 3);

//        delete pipeline_opt;

//        times.saveCSV(outputDir + "/times.txt");

        if(split)
        {
            delete Xtrain;
            delete ytrain;
        }

        return EXIT_SUCCESS;
    }
    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

}


GurlsOptionsList *runPipeline(const gMat2D<T>& Xte, const gMat2D<T>&yte, const gMat2D<T>& Xtr, const gMat2D<T> &alpha, double sigma)
{
    OptTaskSequence *seq = new OptTaskSequence();
    *seq << "predkernel:traintest" << "pred:dual" << "perf:macroavg";

    GurlsOptionsList * process = new GurlsOptionsList("processes", false);

    OptProcess* process1 = new OptProcess();
    *process1 << GURLS::computeNsave << GURLS::computeNsave << GURLS::computeNsave;
    process->addOpt("one", process1);


    GurlsOptionsList* opt = new GurlsOptionsList("ichol", true);
    opt->addOpt("seq", seq);
    opt->addOpt("processes", process);


    GurlsOptionsList *optimizer = new GurlsOptionsList("optimizer");
    GurlsOptionsList* kernel = new GurlsOptionsList("kernel");
    GurlsOptionsList* paramsel = new GurlsOptionsList("paramsel");

    kernel->addOpt("type", "rbf");
    paramsel->addOpt("sigma", new OptNumber(sigma));

//     opt.optimizer.C = alpha_opt:
    gMat2D<T> *alpha_mat = new gMat2D<T>(alpha);
    optimizer->addOpt("C", new OptMatrix<gMat2D<T> >(*alpha_mat));

//     opt.optimizer.X = Xtr:
    gMat2D<T> *Xtr_mat = new gMat2D<T>(Xtr);
    optimizer->addOpt("X", new OptMatrix<gMat2D<T> >(*Xtr_mat));

    opt->addOpt("optimizer", optimizer);
    opt->addOpt("kernel", kernel);
    opt->addOpt("paramsel", paramsel);

    GURLS G;

    std::string jobId0("one");

    G.run(Xte, yte, *opt, jobId0);

    return opt;
}
