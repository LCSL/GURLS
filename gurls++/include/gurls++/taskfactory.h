/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-2013, IIT@MIT Lab
 * All rights reserved.
 *
 * author:  M. Santoro
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

#ifndef _GURLS_TASKFACTORY_H_
#define _GURLS_TASKFACTORY_H_


#include "gurls++/optlist.h"
#include "gurls++/exceptions.h"
#include "gurls++/split.h"
#include "gurls++/confidence.h"
#include "gurls++/paramsel.h"
#include "gurls++/norm.h"
#include "gurls++/kernel.h"
#include "gurls++/predkernel.h"
#include "gurls++/perf.h"
#include "gurls++/optimization.h"
#include "gurls++/pred.h"

namespace gurls
{
	
template<typename T>
class ParamSelectionFactory;
template<typename T>
class OptimizerFactory;
template<typename T>
class SplitFactory;
template<typename T>
class PredictionFactory;
template<typename T>
class PerformanceFactory;
template<typename T>
class KernelFactory;
template<typename T>
class PredKernelFactory;
template<typename T>
class NormFactory;
template<typename T>
class ConfidenceFactory;

template<typename T>
class TaskFactory
{
public:
	TaskFactory(){};
	virtual Task<T>* operator()(const std::string& id) throw(BadTaskCreation) = 0;

	static Task<T> *factory(const std::string &id)
	{
    	if(!initialized())
        	initialize();

    	size_t found = id.find(gurls::TASKDESC_SEPARATOR);

    	if (found==std::string::npos)
        	throw BadTaskCreation(id);

    	std::string category = id.substr(0, found);
    	std::string type = id.substr(found+1);

if(m_registeredFactories().find(category) == m_registeredFactories().end())
        	throw BadTaskCreation(id);


    	typedef TaskFactory<T> FactoryType;
    	typedef std::vector<FactoryType*> FVector;
    	const FVector &factories = m_registeredFactories()[category];

    	for(typename FVector::const_iterator it = factories.begin(),
           end = factories.end(); it != end; ++it)
    	{
        	try{
            	    return (*(*it))(type);
        	}
        	catch(BadTaskCreation){}
    	}

    	throw BadTaskCreation(id);
	}

	static void registerFactory(const std::string &category,
                                   TaskFactory<T> *factory)
	{
    	    m_registeredFactories()[category].push_back(factory);
	}

private:

	static void initialize();

	static bool &initialized()
	{
    	    static bool ret = false;
    	    return ret;
	}

static std::map<std::string, std::vector<TaskFactory<T>*> >  &m_registeredFactories()
	{
    	    static std::map<std::string, std::vector<TaskFactory<T>*> > ret;
    	    return ret;
	}
};

template<typename T>
void TaskFactory<T>::initialize()
{
	m_registeredFactories()["optimizer"].push_back(new OptimizerFactory<T>());
	m_registeredFactories()["paramsel"].push_back(new ParamSelectionFactory<T>());
	m_registeredFactories()["split"].push_back(new SplitFactory<T>());
	m_registeredFactories()["norm"].push_back(new NormFactory<T>());
	m_registeredFactories()["perf"].push_back(new PerformanceFactory<T>());
	m_registeredFactories()["conf"].push_back(new ConfidenceFactory<T>());
	m_registeredFactories()["kernel"].push_back(new KernelFactory<T>());
	m_registeredFactories()["predkernel"].push_back(new PredKernelFactory<T>());
	m_registeredFactories()["pred"].push_back(new PredictionFactory<T>());

	initialized() = true;
}


template<typename T>
class OptimizerFactory: public TaskFactory<T>
{
public:
	OptimizerFactory() :TaskFactory<T>(){}
	gurls::Optimizer<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return Optimizer<T>::factory(id);
	}
};

template<typename T>
class ParamSelectionFactory: public TaskFactory<T>
{
public:
	ParamSelectionFactory() :TaskFactory<T>(){}
	ParamSelection<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return ParamSelection<T>::factory(id);
	}
};

template<typename T>
class SplitFactory: public TaskFactory<T>
{
public:
	SplitFactory() :TaskFactory<T>(){}
	Split<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return Split<T>::factory(id);
	}
};

template<typename T>
class PredictionFactory: public TaskFactory<T>
{
public:
	PredictionFactory() :TaskFactory<T>(){}
	Prediction<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return Prediction<T>::factory(id);
	}
};

template<typename T>
class PerformanceFactory: public TaskFactory<T>
{
public:
	PerformanceFactory() :TaskFactory<T>(){}
	Performance<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return Performance<T>::factory(id);
	}
};

template<typename T>
class PredKernelFactory: public TaskFactory<T>
{
public:
	PredKernelFactory() :TaskFactory<T>(){}
	PredKernel<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return PredKernel<T>::factory(id);
	}
};

template<typename T>
class KernelFactory: public TaskFactory<T>
{
public:
	KernelFactory() :TaskFactory<T>(){}
	Kernel<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return Kernel<T>::factory(id);
	}
};

template<typename T>
class NormFactory: public TaskFactory<T>
{
public:
	NormFactory() :TaskFactory<T>(){}
	Norm<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return Norm<T>::factory(id);
	}
};

template<typename T>
class ConfidenceFactory: public TaskFactory<T>
{
public:
	ConfidenceFactory() :TaskFactory<T>(){}
	Confidence<T>* operator()(const std::string& id) throw(BadTaskCreation)
	{
    	    return Confidence<T>::factory(id);
	}
};
}

#endif // _GURLS_TASKFACTORY_H_

