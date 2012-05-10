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


#include "exceptions.h"
#include "gmat2d.h"
#include "optlist.h"
#include "options.h"

#include "linearkernel.h"

#include "precisionrecall.h"
#include "macroavg.h"

#include "rlsprimal.h"
#include "rlsdual.h"
#include "rlspegasos.h"
#include "loocvprimal.h"
#include "loocvdual.h"
#include "fixlambda.h"
#include "pred.h"
#include "primal.h"
#include "dual.h"

namespace gurls {

class GURLS {

public:
    enum {ignore, compute, computeNsave, load, remove};

    template <typename T>
    void run(const gMat2D<T>& , const gMat2D<T>&,
             GurlsOptionsList&, std::string );

};


template <typename T>
void GURLS::run(const gMat2D<T>& X, const gMat2D<T>& y,
                GurlsOptionsList& opt, std::string processid){

//    time_t start_time,end_time;
    clock_t start_time,end_time;

    Optimizer<T> *taskOpt;
    ParamSelection<T> *taskParSel;
    Prediction<T> *taskPrediction;
    Performance<T> *taskPerformance;
    Kernel<T> *taskKernel;

    try{

        OptTaskSequence* seq = OptTaskSequence::dynacast(opt.getOpt("seq"));
        GurlsOptionsList& processes =
                *GurlsOptionsList::dynacast(opt.getOpt("processes"));

        if (!processes.hasOpt(processid)){
            throw gException(Exception_Gurls_Invalid_ProcessID);
        }
        std::vector<double> process =
                OptNumberList::dynacast( processes.getOpt(processid) )->getValue();

        if (process.size() !=seq->size()){
            throw gException(gurls::Exception_Gurls_Inconsistent_Processes_Number);
        }

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

        if (opt.hasOpt("time")){
            timelist = GurlsOptionsList::dynacast(opt.getOpt("time"));
        }else {
            timelist = new GurlsOptionsList("elapsedtime");
            opt.addOpt("time", timelist);
        }
//        std::vector <double> process_time;
        //process_time.reserve(seq->size());

        std::vector <double> process_time(seq->size(), 0.0);

        //%for i = 1:numel(opt.process) % Go by the length of process.
        //opt.time{jobid} = struct;
        //%end

        std::string reg1;
        std::string reg2;
        std::string fun("");
        std::cout << std::endl
                  <<"####### New task sequence... "
                  << std::endl;

        for (int i = 0; i < seq->size(); i++){

            seq->getTaskAt(i, reg1, reg2);
            std::cout << "\t" << "[Task " << i << ": "
                      << reg1 << "]: " << reg2 << "... ";
            cout.flush();

            switch ( static_cast<int>(process[i]) ) {

            case GURLS::ignore:
//                process_time.assign(i, 0.0);
                std::cout << " ignored." << std::endl;
                break;

            case GURLS::compute:
            case GURLS::computeNsave:
                // WARNING: we should consider the case in which
                // the following statements holds true because the
                // field reg{1} already exists in opt.
                //	case {CPT, CSV, ~isfield(opt,reg{1})}

//                time(&start_time);
                start_time = clock();
                if (!reg1.compare("optimizer")) {
                    taskOpt = Optimizer<T>::factory(reg2);
                    taskOpt->execute(X, y, opt);
                }else if (!reg1.compare("paramsel")) {
                    taskParSel = ParamSelection<T>::factory(reg2);
                    taskParSel->execute(X, y, opt);
                }else if (!reg1.compare("pred")) {
                    taskPrediction = Prediction<T>::factory(reg2);
                    taskPrediction->execute(X, y, opt);
                }else if (!reg1.compare("perf")) {
                    taskPerformance = Performance<T>::factory(reg2);
                    taskPerformance->execute(X, y, opt);
                }else if (!reg1.compare("kernel")) {
                    taskKernel = Kernel<T>::factory(reg2);
                    taskKernel->execute(X, y, opt);
                }
                fun = reg1;
                fun+="_";
                fun+=reg2;
                opt.addOpt(reg1, new OptString(fun));
//                time (&end_time);
                end_time = clock();
                //process_time.assign(i, difftime (end_time,start_time));
                process_time[i] = ((float)end_time - start_time) /((float)CLOCKS_PER_SEC);

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
                std::cout << " skipped." << std::endl;
                break;
            default:
                std::cout << " WARNING: unknown task assignment."
                          << std::endl;
            }

        }

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

        timelist->addOpt(processid, new OptNumberList(process_time));

    }catch (gException& gex) {

        throw gex;
    }

}


}

#include "gurls.hpp"

#endif // _GURLS_GURLS_H_
