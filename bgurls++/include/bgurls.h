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


#ifndef _GURLS_BGURLS_H_
#define _GURLS_BGURLS_H_

#include "gurls.h"
#include "bigarray.h"
#include "bigoptlist.h"

#include "bigparamsel_hoprimal.h"
#include "bigparamsel_calibratesgd.h"

#include "bigoptimizer_rlspegasos.h"
#include "bigoptimizer_rlsprimal.h"

#include "bigperf_macroavg.h"

#include "bigpred_primal.h"

#include "bigsplit_ho.h"

#include <mpi.h>

#include <boost/algorithm/string/erase.hpp>

namespace gurls
{

/**
 * \ingroup Common
 * \brief BGURLS is the class that implements a BGURLS process
 */
class GURLS_EXPORT BGURLS
{

public:

    /**
     * Implements a BGURLS process and stores results of each BGURLS task in opt.
     *
     * \param X input data matrix
     * \param y labels matrix
     * \param opt initial BGURLS options
     * \param processid a job-id number
     *
     */
    template <typename T>
    void run(const BigArray<T>& X, const BigArray<T>& y, GurlsOptionsList& opt, std::string processid, bool hasGurlsProcesses = false);


    template<typename T, class TaskT>
    void runSerialTask(gMat2D<T>& X, gMat2D<T>& y, GurlsOptionsList& opt,
                        int myid, const std::string& taskName, const std::string& reg2, const std::string& dataExchangeFile)
    {
        if(myid == 0)
        {
            TaskT* task = TaskT::factory(reg2);

            GurlsOptionsList* ret = task->execute(X, y, opt);
            ret->save(dataExchangeFile);

            delete task;
            delete ret;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        GurlsOptionsList* ret = new GurlsOptionsList(taskName);
        ret->load(dataExchangeFile);

        opt.removeOpt(taskName);
        opt.addOpt(taskName, ret);

    }

    template<typename T, class TaskT>
    void runBigTask(const BigArray<T>& X, const BigArray<T>& y, GurlsOptionsList& opt,
                    const std::string& taskName, const std::string& reg2)
    {
        TaskT* task = TaskT::factory(reg2);

        GurlsOptionsList* ret = task->execute(X, y, opt);

        delete task;

        opt.removeOpt(taskName);
        opt.addOpt(taskName, ret);
    }

};


template <typename T>
void BGURLS::run(const BigArray<T>& X, const BigArray<T>& y, GurlsOptionsList& opt, std::string processid, bool hasGurlsProcesses)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double begin, end;

//    try{

        OptTaskSequence* seq = OptTaskSequence::dynacast(opt.getOpt("seq"));
        GurlsOptionsList& processes = *GurlsOptionsList::dynacast(opt.getOpt("processes"));

        if (!processes.hasOpt(processid))
            throw gException(Exception_Gurls_Invalid_ProcessID);

        OptProcess* process = processes.getOptAs<OptProcess>(processid);

        if ( process->size() != seq->size())
            throw gException(gurls::Exception_Gurls_Inconsistent_Processes_Number);

        const std::string saveFile = opt.getOptAsString("savefile");
        const std::string sharedDir = opt.getOptAsString("shared_dir");
        const std::string dataExchangeFile = sharedDir + "ret";

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

        GurlsOptionsList *timelist;

        if (opt.hasOpt("time"))
            timelist = GurlsOptionsList::dynacast(opt.getOpt("time"));
        else
        {
            timelist = new GurlsOptionsList("elapsedtime");
            opt.addOpt("time", timelist);
        }


        gMat2D<T>* process_time_vector = new gMat2D<T>(1, seq->size());
        T *process_time = process_time_vector->getData();
        set(process_time, (T)0.0, seq->size());


        std::string reg1;
        std::string reg2;

        if(myid == 0)
            std::cout << std::endl <<"####### New task sequence... " << std::endl;

        gMat2D<T> X_mat;
        gMat2D<T> y_mat;

        if(hasGurlsProcesses)
        {
            X.getMatrix(0, 0, X.rows(), X.cols(), X_mat);
            y.getMatrix(0, 0, y.rows(), y.cols(), y_mat);
        }


        for (unsigned long i = 0; i < seq->size(); ++i)
        {
            MPI_Barrier(MPI_COMM_WORLD);

            seq->getTaskAt(i, reg1, reg2);

            if(myid == 0)
            {
                std::cout << "\t" << "[Task " << i << ": "
                          << reg1 << "]: " << reg2 << "... ";
                std::cout.flush();
            }


            switch( (*process)[i] )
            {
            case GURLS::ignore:
                if(myid == 0)
                    std::cout << " ignored." << std::endl;
                break;

            case GURLS::compute:
            case GURLS::computeNsave:

                begin = MPI_Wtime();

                if (!reg1.compare("optimizer"))
                {
                    runSerialTask<T, Optimizer<T> >(X_mat, y_mat, opt, myid, "optimizer", reg2, dataExchangeFile);
                }
                else if (!reg1.compare("paramsel"))
                {
                    runSerialTask<T, ParamSelection<T> >(X_mat, y_mat, opt, myid, "paramsel", reg2, dataExchangeFile);
                }
                else if (!reg1.compare("pred"))
                {
                    unsigned char isMatrixOption; // old MPI standards doesn't support bool data type

                    if(myid == 0)
                    {
                        Prediction<T> *taskPrediction = Prediction<T>::factory(reg2);

                        GurlsOption* ret = taskPrediction->execute(X_mat, y_mat, opt);

                        isMatrixOption = ret->isA(MatrixOption)? 1 : 0;

                        if(isMatrixOption)
                            OptMatrix<gMat2D<T> >::dynacast(ret)->getValue().save(dataExchangeFile);
                        else
                            GurlsOptionsList::dynacast(ret)->save(dataExchangeFile);

                        delete taskPrediction;
                        delete ret;
                    }

                    MPI_Bcast(&isMatrixOption, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOption* ret;

                    if(isMatrixOption)
                    {
                        gMat2D<T>* mat = new gMat2D<T>();
                        mat->load(dataExchangeFile);

                        ret = new OptMatrix<gMat2D<T> >(*mat);
                    }
                    else
                    {
                        GurlsOptionsList* tmp = new GurlsOptionsList("pred");
                        tmp->load(dataExchangeFile);
                        ret = tmp;
                    }

                    opt.removeOpt("pred");
                    opt.addOpt("pred", ret);
                }
                else if (!reg1.compare("perf"))
                {
                    runSerialTask<T, Performance<T> >(X_mat, y_mat, opt, myid, "perf", reg2, dataExchangeFile);
                }
                else if (!reg1.compare("kernel"))
                {
                    runSerialTask<T, Kernel<T> >(X_mat, y_mat, opt, myid, "kernel", reg2, dataExchangeFile);
                }
                else if (!reg1.compare("norm"))
                {
                    if(myid == 0)
                    {
                        Norm<T> *taskNorm = Norm<T>::factory(reg2);

                        gMat2D<T>* X1 = taskNorm->execute(X_mat, y_mat, opt);

                        delete X1;
                        delete taskNorm;

                        throw gException("Unused return value");

                    }
                }
                else if (!reg1.compare("split"))
                {
                    runSerialTask<T, Split<T> >(X_mat, y_mat, opt, myid, "split", reg2, dataExchangeFile);
                }
                else if (!reg1.compare("predkernel"))
                {
                    runSerialTask<T, PredKernel<T> >(X_mat, y_mat, opt, myid, "predkernel", reg2, dataExchangeFile);
                }
                else if (!reg1.compare("conf"))
                {
                    runSerialTask<T, Confidence<T> >(X_mat, y_mat, opt, myid, "conf", reg2, dataExchangeFile);
                }
                else if (!reg1.compare("bigoptimizer"))
                {
                    runBigTask<T, BigOptimizer<T> >(X, y, opt, "optimizer", reg2);
                }
                else if (!reg1.compare("bigparamsel"))
                {
                    runBigTask<T, BigParamSelection<T> >(X, y, opt, "paramsel", reg2);
                }
                else if (!reg1.compare("bigpred"))
                {
                    runBigTask<T, BigPrediction<T> >(X, y, opt, "pred", reg2);
                }
                else if (!reg1.compare("bigperf"))
                {
                    runBigTask<T, BigPerformance<T> >(X, y, opt, "perf", reg2);
                }
                else if (!reg1.compare("bigsplit"))
                {
                    runBigTask<T, BigSplit<T> >(X, y, opt, "split", reg2);
                }
                else
                    throw gException(Exception_Invalid_TaskSequence);



                end = MPI_Wtime();

                if(myid == 0)
                {
                    process_time[i] = (T)(end-begin);

                    std::cout << " done." << std::endl;
                }

                break;

            case GURLS::load:

                if(loadOpt == NULL)
                    throw gException("Opt savefile not found");

                boost::algorithm::erase_first(reg1, "big");
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

                if(myid ==0)
                    std::cout << " copied" << std::endl;

                break;
            default:
                throw gException("Unknown task assignment");
            }

        }

        if(myid == 0)
        {
            timelist->addOpt(processid, new OptMatrix<gMat2D<T> >(*process_time_vector));

            bool save = false;

            std::cout << std::endl << "Save cycle..." << std::endl;
            for (unsigned long i = 0; i < seq->size(); ++i)
            {
                seq->getTaskAt(i, reg1, reg2);
                boost::algorithm::erase_first(reg1, "big");

                std::cout << "\t" << "[Task " << i << ": " << reg1 << "]: " << reg2 << "... ";
                std::cout.flush();

                switch ( (*process)[i] )
                {
                case GURLS::ignore:
                case GURLS::compute:
                case GURLS::remove:
                    std::cout << "not saved" << std::endl;
                    opt.removeOpt(reg1);
                    break;
                case GURLS::load:
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
        }

        delete loadOpt;

        MPI_Barrier(MPI_COMM_WORLD);

//    }
//    catch (gException& gex)
//    {
//        throw gex;
//    }

}


}

#include "calibratesgd.h"

#endif // _GURLS_BGURLS_H_
