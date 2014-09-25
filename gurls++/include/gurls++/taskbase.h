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

#ifndef _GURLS_TASKBASE_H_
#define _GURLS_TASKBASE_H_

#include "gurls++/exceptions.h"

namespace gurls
{
	
class TaskBase
{
public:
	///
	/// \brief Constructor
	/// \param fieldName The name to be used to index the task outputin the options list
	/// \param taskName The task name
	///
	TaskBase(const std::string &fieldName, const std::string &taskName)
    	:m_fieldName(fieldName), m_taskName(taskName){}
	///
	/// \brief Performs a deep-copy of the task
	///
	virtual TaskBase* clone()
	{
		throw gException(std::string("Undefined 'clone' method in subtask of \"") + fieldName() + std::string("\".") +
			std::string("\nInitialize Task Sequence using standard string, or define 'clone'."));
	}

	///
	/// \brief Gets the name to be used to index the task output in the options list
	///
	const std::string &fieldName()
	{
    	return m_fieldName;
	}

	///
	/// \brief Gets the task name
	///
	const std::string &taskName()
	{
    	return m_taskName;
	}

private:

	/// \brief The name to be used to index the task output in the options list
	std::string m_fieldName;

	/// \brief The task name
	std::string m_taskName;

};

}

#endif // _GURLS_TASKBASE_H_

