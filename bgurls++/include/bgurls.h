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

#include <gurls.h>
#include <bigarray.h>

#include "bigparamsel_hoprimal.h"
#include <mpi/mpi.h>


namespace gurls
{

/**
 * \ingroup Common
 * \brief BGURLS is the class that implements a BGURLS process
 */
class GURLS_EXPORT BGURLS/*: public GURLS*/
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
    void run(const BigArray<T>& X, const BigArray<T>& y, GurlsOptionsList& opt, std::string processid, bool hasGurlsProcesses);

};


template <typename T>
void BGURLS::run(const BigArray<T>& X, const BigArray<T>& y, GurlsOptionsList& opt, std::string processid, bool hasGurlsProcesses)
{

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double begin, end;

    Optimizer<T> *taskOpt;
    ParamSelection<T> *taskParSel;
    Prediction<T> *taskPrediction;
    Performance<T> *taskPerformance;
    Kernel<T> *taskKernel;
    Norm<T> *taskNorm;
    Split<T> *taskSplit;
    PredKernel<T> *taskPredKernel;
    Confidence<T> *taskConfidence;

    BigParamSelection<T> *taskBigParSel;

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
        {
            std::cout << std::endl
                      <<"####### New task sequence... "
                      << std::endl;
        }

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
                    if(myid == 0)
                    {
                        taskOpt = Optimizer<T>::factory(reg2);

                        GurlsOptionsList* ret = taskOpt->execute(X_mat, y_mat, opt);
                        ret->save(dataExchangeFile);

                        delete taskOpt;
                        delete ret;
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOptionsList* ret = new GurlsOptionsList("optimizer");
                    ret->load(dataExchangeFile);

                    opt.removeOpt("optimizer");
                    opt.addOpt("optimizer", ret);
                }
                else if (!reg1.compare("paramsel"))
                {
                    if(myid == 0)
                    {
                        taskParSel = ParamSelection<T>::factory(reg2);

                        GurlsOptionsList* ret = taskParSel->execute(X_mat, y_mat, opt);
                        ret->save(dataExchangeFile);

                        delete taskParSel;
                        delete ret;
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOptionsList* ret = new GurlsOptionsList("paramsel");
                    ret->load(dataExchangeFile);

                    opt.removeOpt("paramsel");
                    opt.addOpt("paramsel", ret);
                }
                else if (!reg1.compare("pred"))
                {
                    unsigned char isMatrixOption; // old MPI standards doesn't support bool data type

                    if(myid == 0)
                    {
                        taskPrediction = Prediction<T>::factory(reg2);

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
                        GurlsOptionsList* ret = new GurlsOptionsList("pred");
                        ret->load(dataExchangeFile);
                    }

                    opt.removeOpt("pred");
                    opt.addOpt("pred", ret);
                }
                else if (!reg1.compare("perf"))
                {
                    if(myid == 0)
                    {
                        taskPerformance = Performance<T>::factory(reg2);

                        GurlsOptionsList* ret = taskPerformance->execute(X_mat, y_mat, opt);
                        ret->save(dataExchangeFile);

                        delete taskPerformance;
                        delete ret;
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOptionsList* ret = new GurlsOptionsList("perf");
                    ret->load(dataExchangeFile);

                    opt.removeOpt("perf");
                    opt.addOpt("perf", ret);

                }
                else if (!reg1.compare("kernel"))
                {
                    if(myid == 0)
                    {
                        taskKernel = Kernel<T>::factory(reg2);

                        GurlsOptionsList* ret = taskKernel->execute(X_mat, y_mat, opt);
                        ret->save(dataExchangeFile);

                        delete taskKernel;
                        delete ret;
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOptionsList* ret = new GurlsOptionsList("kernel");
                    ret->load(dataExchangeFile);

                    opt.removeOpt("kernel");
                    opt.addOpt("kernel", ret);
                }
                else if (!reg1.compare("norm"))
                {
                    if(myid == 0)
                    {
                        taskNorm = Norm<T>::factory(reg2);

                        gMat2D<T>* X1 = taskNorm->execute(X_mat, y_mat, opt);

                        delete X1;
                        delete taskNorm;

                        throw gException("Unused return value");

                    }
                }
                else if (!reg1.compare("split"))
                {
                    if(myid == 0)
                    {
                        taskSplit = Split<T>::factory(reg2);

                        GurlsOptionsList* ret = taskSplit->execute(X_mat, y_mat, opt);
                        ret->save(dataExchangeFile);

                        delete taskSplit;
                        delete ret;
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOptionsList* ret = new GurlsOptionsList("split");
                    ret->load(dataExchangeFile);

                    opt.removeOpt("split");
                    opt.addOpt("split", ret);
                }
                else if (!reg1.compare("predkernel"))
                {
                    if(myid == 0)
                    {
                        taskPredKernel = PredKernel<T>::factory(reg2);

                        GurlsOptionsList* ret = taskPredKernel->execute(X_mat, y_mat, opt);
                        ret->save(dataExchangeFile);

                        delete taskPredKernel;
                        delete ret;
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOptionsList* ret = new GurlsOptionsList("predkernel");
                    ret->load(dataExchangeFile);

                    opt.removeOpt("predkernel");
                    opt.addOpt("predkernel", ret);
                }
                else if (!reg1.compare("conf"))
                {
                    if(myid == 0)
                    {
                        taskConfidence = Confidence<T>::factory(reg2);

                        GurlsOptionsList* ret = taskConfidence->execute(X_mat, y_mat, opt);
                        ret->save(dataExchangeFile);

                        delete taskConfidence;
                        delete ret;
                    }

                    MPI_Barrier(MPI_COMM_WORLD);

                    GurlsOptionsList* ret = new GurlsOptionsList("conf");
                    ret->load(dataExchangeFile);

                    opt.removeOpt("conf");
                    opt.addOpt("conf", ret);

                }
                else if (!reg1.compare("bigparamsel"))
                {
                    taskBigParSel = BigParamSelection<T>::factory(reg2);

                    GurlsOptionsList* ret = taskBigParSel->execute(X, y, opt);

                    opt.removeOpt("paramsel");
                    opt.addOpt("paramsel", ret);

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

                cout << "Load " << myid << "loadopt=" << loadOpt;
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

        if(myid == 0)
        {
            timelist->addOpt(processid, new OptMatrix<gMat2D<T> >(*process_time_vector));

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
