/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, Matteo Santoro
 * All rights reserved.
 *
 * author:  M. Santoro
 * email:   matteo.santoro@gmail.com
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


#ifndef _GURLS_GURLS_H_
#define _GURLS_GURLS_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <exception>
#include <ctime>

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include "exports.h"
#include "exceptions.h"
#include "gmat2d.h"
#include "optlist.h"
#include "options.h"
#include "optarray.h"
#include "optfunction.h"
#include "optmatrix.h"

#include "linearkernel.h"
#include "rbfkernel.h"
#include "chisquaredkernel.h"

#include "predkerneltraintest.h"

#include "precisionrecall.h"
#include "macroavg.h"
#include "rmse.h"

#include "rlsauto.h"
#include "rlsprimal.h"
#include "rlsprimalr.h"
#include "rlsdual.h"
#include "rlsdualr.h"
#include "rlspegasos.h"
#include "rlsgp.h"

#include "loocvprimal.h"
#include "loocvdual.h"
#include "fixlambda.h"
#include "fixsiglam.h"
#include "siglam.h"
#include "siglamho.h"
#include "hoprimal.h"
#include "hodual.h"

#include "hogpregr.h"
#include "loogpregr.h"
#include "siglamhogpregr.h"
#include "siglamloogpregr.h"

#include "pred.h"
#include "primal.h"
#include "dual.h"
#include "predgp.h"

#include "norml2.h"
#include "normtestzscore.h"
#include "normzscore.h"

#include "splitho.h"

#include "boltzman.h"
#include "boltzmangap.h"
#include "gap.h"
#include "maxscore.h"

namespace gurls {

    /**
     * \ingroup Common
     * \brief GURLS is the class that implements a GURLS process
     */

    class GURLS_EXPORT GURLS
    {

    public:

        /**
         * Execution options for a GURLS task
         */
//        enum Action {ignore, compute, computeNsave, load, remove};

        typedef OptProcess::Action Action;

        static const Action ignore = OptProcess::ignore;
        static const Action compute = OptProcess::compute;
        static const Action computeNsave = OptProcess::computeNsave;
        static const Action load = OptProcess::load;
        static const Action remove = OptProcess::remove;

        /**
         * Implements a GURLS process and stores results of each GULRS task in opt.
         *
         * \param X input data matrix
         * \param y labels matrix
         * \param opt initial GURLS options
         * \param processid a job-id number
         *
         */
        template <typename T>
        void run(const gMat2D<T>& X, const gMat2D<T>& y,
                 GurlsOptionsList& opt, std::string processid);

};


template <typename T>
void GURLS::run(const gMat2D<T>& X, const gMat2D<T>& y,
                GurlsOptionsList& opt, std::string processid)
{

    boost::posix_time::ptime begin, end;
    boost::posix_time::time_duration diff;

    Optimizer<T> *taskOpt;
    ParamSelection<T> *taskParSel;
    Prediction<T> *taskPrediction;
    Performance<T> *taskPerformance;
    Kernel<T> *taskKernel;
    Norm<T> *taskNorm;
    Split<T> *taskSplit;
    PredKernel<T> *taskPredKernel;
    Confidence<T> *taskConfidence;

//    try{

        OptTaskSequence* seq = OptTaskSequence::dynacast(opt.getOpt("seq"));
        GurlsOptionsList& processes = *GurlsOptionsList::dynacast(opt.getOpt("processes"));

        if (!processes.hasOpt(processid))
            throw gException(Exception_Gurls_Invalid_ProcessID);

//        std::vector<double> process = OptNumberList::dynacast( processes.getOpt(processid) )->getValue();
        OptProcess* process = processes.getOptAs<OptProcess>(processid);

//        if ((long)process.size() != seq->size())
        if ( process->size() != seq->size())
            throw gException(gurls::Exception_Gurls_Inconsistent_Processes_Number);

        const std::string saveFile = opt.getOptAsString("savefile");

        GurlsOptionsList* loadOpt = new GurlsOptionsList("load");
        try
        {
            loadOpt->load(saveFile);
        }
        catch(gException & ex)
        {
            delete loadOpt;
            loadOpt = NULL;
        }

        GurlsOption* tmpOpt;

        //% Load and copy

        //if exist(opt.savefile) == 2
        //	t = load(opt.savefile);
        //	if isfield(t.opt,'time');
        //		opt.time = t.opt.time;
        //	end
        //else
        //	fprintf('Could not load %s. Starting from scratch.\n', opt.savefile);
        //end
        //%try
        //%	t = load(opt.savefile);
        //%	if isfield(t.opt,'time')
        //%		opt.time = t.opt.time;
        //%	end
        //%catch
        //%	fprintf('Could not load %s. Starting from scratch.\n', opt.savefile);
        //%end

        GurlsOptionsList *timelist;

        if (opt.hasOpt("time"))
            timelist = GurlsOptionsList::dynacast(opt.getOpt("time"));
        else
        {
            timelist = new GurlsOptionsList("elapsedtime");
            opt.addOpt("time", timelist);
        }


//        std::vector <double> process_time(seq->size(), 0.0);
        gMat2D<T>* process_time_vector = new gMat2D<T>(1, seq->size());
        T *process_time = process_time_vector->getData();
        set(process_time, (T)0.0, seq->size());

        //%for i = 1:numel(opt.process) % Go by the length of process.
        //opt.time{jobid} = struct;
        //%end

        std::string reg1;
        std::string reg2;
//        std::string fun("");
        std::cout << std::endl
                  <<"####### New task sequence... "
                  << std::endl;

        for (unsigned long i = 0; i < seq->size(); ++i)
        {
            seq->getTaskAt(i, reg1, reg2);

            std::cout << "\t" << "[Task " << i << ": "
                      << reg1 << "]: " << reg2 << "... ";
            std::cout.flush();


//            switch ( static_cast<int>(process[i]) )
            switch( (*process)[i] )
            {
            case GURLS::ignore:
                std::cout << " ignored." << std::endl;
                break;

            case GURLS::compute:
            case GURLS::computeNsave:
                // WARNING: we should consider the case in which
                // the following statements holds true because the
                // field reg{1} already exists in opt.
                //	case {CPT, CSV, ~isfield(opt,reg{1})}

                begin = boost::posix_time::microsec_clock::local_time();

                if (!reg1.compare("optimizer"))
                {
                    taskOpt = Optimizer<T>::factory(reg2);
                    GurlsOption* ret = taskOpt->execute(X, y, opt);
                    opt.removeOpt("optimizer");
                    opt.addOpt("optimizer", ret);
                }
                else if (!reg1.compare("paramsel"))
                {
                    taskParSel = ParamSelection<T>::factory(reg2);
                    GurlsOption* ret = taskParSel->execute(X, y, opt);
                    opt.removeOpt("paramsel");
                    opt.addOpt("paramsel", ret);
                }
                else if (!reg1.compare("pred"))
                {
                    taskPrediction = Prediction<T>::factory(reg2);
                    GurlsOption* ret = taskPrediction->execute(X, y, opt);
                    opt.removeOpt("pred");
                    opt.addOpt("pred", ret);
                }
                else if (!reg1.compare("perf"))
                {
                    taskPerformance = Performance<T>::factory(reg2);
                    GurlsOption* ret = taskPerformance->execute(X, y, opt);
                    opt.removeOpt("perf");
                    opt.addOpt("perf", ret);
                }
                else if (!reg1.compare("kernel"))
                {
                    taskKernel = Kernel<T>::factory(reg2);
                    GurlsOption* ret = taskKernel->execute(X, y, opt);
                    opt.removeOpt("kernel");
                    opt.addOpt("kernel", ret);
                }
                else if (!reg1.compare("norm"))
                {
                    taskNorm = Norm<T>::factory(reg2);
                    gMat2D<T>* X1 = taskNorm->execute(X, y, opt);
                    delete X1;
                    throw gException("Unused return value");
                }
                else if (!reg1.compare("split"))
                {
                    taskSplit = Split<T>::factory(reg2);
                    GurlsOption* ret = taskSplit->execute(X, y, opt);
                    opt.removeOpt("split");
                    opt.addOpt("split", ret);
                }
                else if (!reg1.compare("predkernel"))
                {
                    taskPredKernel = PredKernel<T>::factory(reg2);
                    GurlsOption* ret = taskPredKernel->execute(X, y, opt);
                    opt.removeOpt("predkernel");
                    opt.addOpt("predkernel", ret);
                }
                else if (!reg1.compare("conf"))
                {
                    taskConfidence = Confidence<T>::factory(reg2);
                    GurlsOption* ret = taskConfidence->execute(X, y, opt);
                    opt.removeOpt("conf");
                    opt.addOpt("conf", ret);
                }

//                fun = reg1;
//                fun+="_";
//                fun+=reg2;
//                opt.addOpt(reg1, new OptString(fun));

                end = boost::posix_time::microsec_clock::local_time();
                diff = end-begin;

                process_time[i] = ((T)diff.total_milliseconds())/1000.0;

                //		fName = [reg{1} '_' reg{2}];
                //		fun = str2func(fName);
                //		tic;
                //		opt = setfield(opt, reg{1}, fun(X, y, opt));
                //		opt.time{jobid} = setfield(opt.time{jobid},reg{1}, toc);
                //		fprintf('\tdone\n');
                std::cout << " done." << std::endl;
                break;

            case GURLS::load:
                //	case LDF,
                //		if exist('t','var') && isfield (t.opt, reg{1})
                //			opt = setfield(opt, reg{1}, getfield(t.opt, reg{1}));
                //			fprintf('\tcopied\n');
                //		else
                //			fprintf('\tcopy failed\n');
                //		end
//                std::cout << " skipped." << std::endl;

                if(loadOpt == NULL)
                    throw gException("Opt savefile not found");
                if(!loadOpt->hasOpt(reg1))
                {
                    std::string s = "Task " + reg1 + " not found in opt savefile";
                    gException e(s);
                    throw e;
                }

                opt.removeOpt(reg1);
                tmpOpt = loadOpt->getOpt(reg1);
                loadOpt->removeOpt(reg1, false);
                opt.addOpt(reg1, tmpOpt);
                std::cout << " copied" << std::endl;

                break;
            default:
                throw gException("Unknown task assignment");
            }

        }

//        timelist->addOpt(processid, new OptNumberList(process_time));
        timelist->addOpt(processid, new OptMatrix<gMat2D<T> >(*process_time_vector));

        //fprintf('\nSave cycle...\n');
        //% Delete whats not necessary
        //for i = 1:numel(process)
        //	fprintf('[Job %d: %15s] %15s: ',jobid, reg{1}, reg{2});
        //	reg = regexp(seq{i},':','split');
        //	switch process(i)
        //		case {CSV, LDF}
        //			fprintf('\tsaving..\n');
        //		otherwise
        //			if isfield (opt, reg{1})
        //				opt = rmfield(opt, reg{1});
        //				fprintf('\tremoving..\n');
        //			else
        //				fprintf('\tnot found..\n');
        //			end
        //	end
        //end
        //save(opt.savefile, 'opt', '-v7.3');
        //fprintf('Saving opt in %s\n', opt.savefile);

        bool save = false;

        std::cout << std::endl << "Save cycle..." << std::endl;
        for (unsigned long i = 0; i < seq->size(); ++i)
        {
            seq->getTaskAt(i, reg1, reg2);
            std::cout << "\t" << "[Task " << i << ": " << reg1 << "]: " << reg2 << "... ";
            std::cout.flush();

            switch ( (*process)[i] )
            {
            case GURLS::ignore:
            case GURLS::compute:
            case GURLS::load:
            case GURLS::remove:
                std::cout << "not saved" << std::endl;
                opt.removeOpt(reg1);
                break;
            case GURLS::computeNsave:
                std::cout << " saving" << std::endl;
                save = true;
                break;
            }
        }

        if(save)
        {
            std::cout << std::endl << "Saving opt in " << saveFile << std::endl;
            opt.save(saveFile);
        }

        delete loadOpt;

//    }
//    catch (gException& gex)
//    {
//        throw gex;
//    }

}


}

#include "gurls.hpp"
#include "calibratesgd.h"

#endif // _GURLS_GURLS_H_
