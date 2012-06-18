#ifndef _GURLS_TESTUTILS_H_
#define _GURLS_TESTUTILS_H_

#include <fstream>
#include <string>
#include <cassert>

#include "exceptions.h"
#include "gmat2d.h"

#include "optlist.h"
#include "gurls.h"

using namespace gurls;

void printHelp(int argc, char*argv[])
{
    std::cout << "========================================================================================"<< std::endl
    << " Wrong parameters number ("<<argc <<")." << std::endl
    << " Provide a valid path for training, test and output files using the following syntax:" << std::endl
    << " \n\n\t " << argv[0] << " xtr xte ytr yte times performances [saveFile]" << std::endl
    << "========================================================================================" << std::endl << std::endl;
}

template<typename T>
GurlsOptionsList* runTest(std::string& name, OptTaskSequence &seq, std::vector<OptNumberList*> &processes, std::vector<std::string> &jobIds,
             const gMat2D<T>& Xtr, const gMat2D<T>& ytr, const gMat2D<T>& Xte, const gMat2D<T>& yte)
{
    GurlsOptionsList* opt = new GurlsOptionsList(name, true);

    opt->addOpt("seq", &seq);
    opt->addOpt("hoperf", "macroavg");

    GurlsOptionsList *process = new GurlsOptionsList("processes", false);

    assert(processes.size() == jobIds.size());

    for(unsigned int i=0; i< jobIds.size(); ++i)
        process->addOpt(jobIds[i], processes[i]);

    opt->addOpt("processes", process);


    if(!processes.empty())
    {
        GURLS G;

        G.run(Xtr, ytr, *opt, jobIds[0]);

        for(unsigned int i = 1; i< jobIds.size(); ++i)
            G.run(Xte, yte, *opt, jobIds[i]);
    }

    return opt;
}

void saveTimes( const std::string &fileName,
                GurlsOptionsList &opt,
                const std::string &tasks,
                const std::vector<OptNumberList*> &processes,
                const std::vector<std::string> &jobIds)
{
    std::ofstream str(fileName.c_str());

    if(!str.is_open())
        throw gException("Could not open file " + fileName);

    GurlsOptionsList *timelist = GurlsOptionsList::dynacast(opt.getOpt("time"));
    std::vector<std::vector<double>* > times;

    for(unsigned int i=0; i<jobIds.size(); ++i)
        times.push_back( &(OptNumberList::dynacast(timelist->getOpt(jobIds[i]))->getValue()));


    str << tasks << std::endl;

    for(unsigned int i=1; i<jobIds.size(); ++i)
    {
        std::vector<double>::const_iterator p0_it = processes[0]->getValue().begin();
        std::vector<double>::const_iterator p0_end = processes[0]->getValue().end();

        std::vector<double>::const_iterator t0_it = times[0]->begin();
        std::vector<double>::const_iterator ti_it = times[i]->begin();


        for(; p0_it != p0_end; ++p0_it, ++t0_it, ++ti_it)
        {
            if(static_cast<int>(*p0_it) == GURLS::computeNsave)
                str << *t0_it << " ";
            else
                str << *ti_it << " ";
        }
        str << std::endl;
    }

    str.close();
}

template<typename T>
void savePerformances(const std::string &fileName, GurlsOptionsList &opt)
{
    std::ofstream str(fileName.c_str());

    if(!str.is_open())
        throw gException("Could not open file " + fileName);

    GurlsOptionsList *perf = static_cast<GurlsOptionsList*>(opt.getOpt("perf"));
    GurlsOption *acc_opt = perf->getOpt("acc");
    gMat2D<T> &acc_mat = (OptMatrix<gMat2D<T> >::dynacast(acc_opt))->getValue();

    for(T *it = acc_mat.getData(), *end = acc_mat.getData()+acc_mat.getSize(); it != end; ++it)
        str << *it << " ";

    str << std::endl;
    str.close();
}

template<typename T>
gMat2D<T> * readFile(const std::string &fileName, bool ROWM = true )
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

#endif //_GURLS_TESTUTILS_H_
