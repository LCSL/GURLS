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

#ifndef _GURLS_OPTTASK_H_
#define _GURLS_OPTTASK_H_

#include <gurls++/task.h>
#include <gurls++/taskfactory.h>

namespace gurls
{
	
class OptTask :public GurlsOption
{
public:
	///
  	/// Constructor from string
  	///
	OptTask(const std::string &name): GurlsOption(TaskOption), m_task(NULL), m_taskString(name){}

	///
  	/// Constructor from TaskBase
  	///
	OptTask(TaskBase* task):GurlsOption(TaskOption), m_task(task)
	{
    	m_taskString = std::string(task->fieldName() + ":" + task->taskName());
	}

	///
  	/// Copy constructor
  	///
	OptTask(const OptTask& task): GurlsOption(TaskOption), m_task(NULL), m_taskString(task.m_taskString)
	{
		if (task.m_task != NULL)
			m_task = task.m_task->clone();
	}

	///
  	/// Destructor
  	///
	~OptTask()
	{
    	if(m_task != NULL)
        	delete m_task;
	}

	///
  	/// Returns an umanaged Task<T>*
  	///
	template<typename T>
	Task<T> *getValue() const
	{
    	if(m_task == NULL)
    		return dynamic_cast<Task<T>*>(TaskFactory<T>::factory(m_taskString));
    	return dynamic_cast<Task<T>*>(m_task->clone());
	}

	///
	/// Returns an unmanaged Task<T>* containing a task copied from opt.optionName
	///
	template<typename T, class TaskType>
	static Task<T> *getValue(const GurlsOptionsList &opt, const std::string &optionName)
	{
    	const gurls::GurlsOption *option = opt.getOpt(optionName);
    	return option->isA(TaskOption)?
                	OptTask::dynacast(option)->getValue<T>():
                	TaskType::factory(OptString::dynacast(option)->getValue());
	}

	///
  	/// Checks if the option has the given type
  	///
	virtual bool isA(OptTypes id) const
	{return id==TaskOption;}

	///
  	/// Tries to cast a pointer to a generic option to a pointer to
	/// an \ref OptTask
  	///
	static OptTask* dynacast(GurlsOption* opt)
	{
	if (opt->isA(TaskOption))
		return static_cast<OptTask*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
	}

	///
  	/// Tries to cast a pointer to a generic option to a pointer to
	/// an \ref OptTask
  	///
	static const OptTask* dynacast(const GurlsOption* opt)
	{
    if (opt->isA(TaskOption))
        return static_cast<const OptTask*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
	}
	
	///
  	/// Returns a const reference to m_taskString
  	///
	const std::string& getString(){return m_taskString;}
	
	///
  	/// Returns a const pointer to m_task
  	///
	const TaskBase* getTask(){return m_task;}
	
	///
  	/// Writes an OptTask to stream
  	///
	virtual std::ostream& operator<<(std::ostream& os) const
	{
    	os << m_taskString;
    	return os;
	}

protected:
	TaskBase *m_task;
	std::string m_taskString;
};
}

#endif // _GURLS_OPTTASK_H_

