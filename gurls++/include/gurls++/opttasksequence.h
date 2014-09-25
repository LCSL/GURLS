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

#ifndef _GURLS_OPTTASKSEQ_H_
#define _GURLS_OPTTASKSEQ_H_

#include "gurls++/taskbase.h"
#include "gurls++/task.h"
#include "gurls++/opttask.h"
#include "gurls++/optlist.h"
#include "gurls++/exceptions.h"

namespace gurls
{
	/**
  * \ingroup Settings
  * \brief OptTaskSequence is an option containing
  * a sequence of task that forms a pipeline
  */

class OptTask;
class GURLS_EXPORT OptTaskSequence: public GurlsOption
{
public:

	typedef std::vector<OptTask*> ValueType;

	///
	/// Empty constructor
	///
	OptTaskSequence():GurlsOption(TaskSequenceOption){}

	///
	///Destructor
	///
	~OptTaskSequence()
	{
    	typedef ValueType::const_iterator TaskIterator;
    	for(TaskIterator it = m_tasks.begin(), end = m_tasks.end(); it != end; ++it)
        	delete *it;
	}
	
	///
	///Copy constructor
	///
	OptTaskSequence(const OptTaskSequence &other):GurlsOption(TaskSequenceOption){
		for(unsigned long i=0; i<other.size(); ++i)
			(*this)<<new OptTask(*(other.getValue()[i]));
	}

	///
	/// Copies the matrix from an existing \ref OptTaskSequence
	///
	OptTaskSequence& operator=(const OptTaskSequence& other)
	{
		for(unsigned int i=0; i<other.size(); ++i)
			(*this)<<new OptTask(*(other.getValue()[i]));
	}

	///
	/// Checks if the option has the given type
	///
	virtual bool isA(OptTypes id) const;
	
	///
	///Checks if a string is in the right format to be passed to an OptTaskSequence, copies the two substrings into type and name
	///
	static bool isValid(const std::string & str, std::string& type, std::string& name);

	///
	/// Tries to cast a pointer to a generic option to a pointer to an \ref OptTaskSequence
	///
	static OptTaskSequence* dynacast(GurlsOption* opt);

	///
	/// Tries to cast a pointer to a generic option to a pointer to an \ref OptTaskSequence
	///
	static const OptTaskSequence* dynacast(const GurlsOption* opt);

	///
	/// Gets the task at a given index
	///
	OptTask* operator[](int index) const
	{
		return m_tasks[index];
	}
	
	///
	/// Gets the task at a given index, in string format
	///
	void getTaskAt(int index, std::string &taskdesc, std::string &taskname)
	{
    if (!isValid(m_tasks[index]->getString(), taskdesc, taskname))
        throw gException(Exception_Invalid_TaskSequence);
	}

	///
	/// Adds a task to the list
	///
	void addTask(const std::string newtask);

	///
	/// Adds a task to the list
	///
	OptTaskSequence& operator<<(const std::string& str)
	{
    	m_tasks.push_back(new OptTask(str));
    	return *this;
	}

	///
	/// Adds a task to the list
	///
	OptTaskSequence& operator<<(OptTask *task)
	{
    	m_tasks.push_back(task);
    	return *this;
	}

	///
	/// Clears the sequence
	///
	void clear();

	///
	/// Adds a task to the list
	///
	OptTaskSequence& operator<<(TaskBase *task)
	{
    	m_tasks.push_back(new OptTask(task));
    	return *this;
	}

	///
	/// Returns the entire vector
	///
	ValueType& getValue()
	{
		return m_tasks;
	}
	
	///
	/// Returns the entire vector
	///
    const ValueType& getValue() const
    {
        return m_tasks;
    }

	///
	/// Returns the number of tasks into the sequence
	///
	unsigned long size() const;
	
	///
	/// Writes the OptTaskSequence to stream
	///
	virtual std::ostream& operator<<(std::ostream& os) const
	{
		os<<"Sequence of size:"<< size();
    	return os;
	}

	class iterator
	{
	public:
    	iterator(ValueType::iterator it):m_it(it){}

    	inline void operator =(iterator &other){ m_it = other.m_it;}
    	inline bool operator !=(iterator &other) const {return m_it != other.m_it;}
    	inline void operator++(){++m_it;}

    	OptTask &operator*()
    	{
        	return **m_it;
    	}

	private:
    	ValueType::iterator m_it;
	};

	///
	/// Returns an iterator to the beginning of the sequence
	///
	iterator begin(){return iterator(m_tasks.begin());}

	///
	/// Returns an iterator to the end of the sequence
	///
	iterator end(){return iterator(m_tasks.end());}

protected:
	ValueType m_tasks;	
};
}

#endif // _GURLS_OPTTASKSEQ_H_

