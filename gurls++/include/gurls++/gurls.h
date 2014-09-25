/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-1013, Matteo Santoro
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

#include "gurls++/exports.h"
#include "gurls++/exceptions.h"
#include "gurls++/gmat2d.h"
#include "gurls++/optlist.h"
#include "gurls++/options.h"
#include "gurls++/optarray.h"
#include "gurls++/optfunction.h"
#include "gurls++/optmatrix.h"
#include "gurls++/opttask.h"
#include "gurls++/opttasksequence.h"

#include "gurls++/linearkernel.h"
#include "gurls++/rbfkernel.h"
#include "gurls++/chisquaredkernel.h"

#include "gurls++/predkerneltraintest.h"

#include "gurls++/precisionrecall.h"
#include "gurls++/macroavg.h"
#include "gurls++/totalavg.h"
#include "gurls++/rmse.h"
#include "gurls++/abserr.h"

#include "gurls++/rlsauto.h"
#include "gurls++/rlsprimal.h"
#include "gurls++/rlsprimalr.h"
#include "gurls++/rlsdual.h"
#include "gurls++/rlsdualr.h"
#include "gurls++/rlspegasos.h"
#include "gurls++/rlsgp.h"
#include "gurls++/rlsprimalrecinit.h"
#include "gurls++/rlsprimalrecupdate.h"
#include "gurls++/rlsrandfeats.h"
#include "gurls++/rlsprimalrecinitcholesky.h"
#include "gurls++/rlsprimalrecupdatecholesky.h"

#include "gurls++/loocvprimal.h"
#include "gurls++/loocvdual.h"
#include "gurls++/fixlambda.h"
#include "gurls++/fixsiglam.h"
#include "gurls++/siglam.h"
#include "gurls++/siglamho.h"
#include "gurls++/hoprimal.h"
#include "gurls++/hodual.h"

#include "gurls++/hogpregr.h"
#include "gurls++/loogpregr.h"
#include "gurls++/siglamhogpregr.h"
#include "gurls++/siglamloogpregr.h"

#include "gurls++/pred.h"
#include "gurls++/primal.h"
#include "gurls++/dual.h"
#include "gurls++/predgp.h"
#include "gurls++/predrandfeats.h"

#include "gurls++/norml2.h"
#include "gurls++/normtestzscore.h"
#include "gurls++/normzscore.h"

#include "gurls++/splitho.h"

#include "gurls++/boltzman.h"
#include "gurls++/boltzmangap.h"
#include "gurls++/gap.h"
#include "gurls++/maxscore.h"
 
#include "gurls++/wrapper.h"

#include "gurls++/taskfactory.h"

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

    boost::posix_time::ptime begintime, endtime;
    boost::posix_time::time_duration diff;


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
        catch(gException & /*ex*/)
        {
            delete loadOpt;
            loadOpt = NULL;
        }

        GurlsOption* tmpOpt;

        GurlsOptionsList *timelist;
		
		std::string fieldName;
		std::string taskName;

        if (opt.hasOpt("time"))
            timelist = GurlsOptionsList::dynacast(opt.getOpt("time"));
        else
        {
            timelist = new GurlsOptionsList("elapsedtime");
            opt.addOpt("time", timelist);
        }

		int todisk=1;
        if (opt.hasOpt("todisk"))
			todisk=(int)opt.getOptAsNumber("todisk");

//        std::vector <double> process_time(seq->size(), 0.0);
        gMat2D<T>* process_time_vector = new gMat2D<T>(1, seq->size());
        T *process_time = process_time_vector->getData();
        set(process_time, (T)0.0, seq->size());

        //%for i = 1:numel(opt.process) % Go by the length of process.
        //opt.time{jobid} = struct;
        //%end

//        std::string fun("");
        std::cout << std::endl
                  <<"####### New task sequence... "
                  << std::endl;
		
	unsigned int i=0;
	for(gurls::OptTaskSequence::iterator it = seq->begin(), end = seq->end(); it != end; ++it, ++i)
        {
        gurls::OptTask &taskOption = *it;
		gurls::Task<T>* task = taskOption.getValue<T>();
		OptTaskSequence::isValid(taskOption.getString(), fieldName, taskName);
		
    	std::cout << "\t" << "[Task " << i << ": "
              	<< fieldName << "]: " << taskName << "... ";
    	std::cout.flush();



//            switch ( static_cast<int>(process[i]) )
            switch( (*process)[i] )
            {
            case GURLS::ignore:
                std::cout << " ignored." << std::endl;
                break;

            case GURLS::compute:
            case GURLS::computeNsave:
                begintime = boost::posix_time::microsec_clock::local_time();
             
    			opt.removeOpt(task->fieldName());
    			opt.addOpt(task->fieldName(), task->execute(X, y, opt));

                endtime = boost::posix_time::microsec_clock::local_time();
                diff = endtime-begintime;

                process_time[i] = ((T)diff.total_milliseconds())/1000.0;

                std::cout << " done." << std::endl;
                break;

            case GURLS::load:

                if(loadOpt == NULL)
				{
					delete task;
                    throw gException("Opt savefile not found");
				}
                if(!loadOpt->hasOpt(task->fieldName()))
                {
                    std::string s = "Task " +  task->fieldName() + " not found in opt savefile";
                    gException e(s);
					delete task;
                    throw e;
                }

                opt.removeOpt( task->fieldName());
                tmpOpt = loadOpt->getOpt( task->fieldName());
                loadOpt->removeOpt( task->fieldName(), false);
                opt.addOpt( task->fieldName(), tmpOpt);
                std::cout << " copied" << std::endl;

                break;
            default:
				delete task;
                throw gException("Unknown task assignment");
            }
			
			delete task;
        }

        timelist->removeOpt(processid);
        timelist->addOpt(processid, new OptMatrix<gMat2D<T> >(*process_time_vector));

        bool save = false;

        std::cout << std::endl << "Save cycle..." << std::endl;
	i=0;
	for(gurls::OptTaskSequence::iterator it = seq->begin(),
												end = seq->end(); it != end; ++it, ++i)
        {
			gurls::OptTask &taskOption = *it;

			OptTaskSequence::isValid(taskOption.getString(), fieldName, taskName);
		
    		std::cout << "\t" << "[Task " << i << ": "<< fieldName << "]: " << taskName << "... ";
    		std::cout.flush();

            switch ( (*process)[i] )
            {
            case GURLS::ignore:
            case GURLS::compute:
            case GURLS::remove:
                std::cout << "not saved" << std::endl;
                opt.removeOpt(fieldName);
                break;
            case GURLS::load:
            case GURLS::computeNsave:
                std::cout << " saving" << std::endl;
                save = true;
                break;
            }
        }

        if(save && todisk)
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
