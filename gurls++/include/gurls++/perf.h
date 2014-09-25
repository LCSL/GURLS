/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
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


#ifndef _GURLS_PERF_H_
#define _GURLS_PERF_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "gurls++/gmath.h"
#include "gurls++/options.h"
#include "gurls++/optlist.h"
#include "gurls++/exceptions.h"


namespace gurls
{
	
template <typename T>
class PerfMacroAvg;

template <typename T>
class PerfTotalAvg;

template <typename T>
class PerfPrecRec;

template <typename T>
class PerfRmse;

template <typename T>
class PerfAbsErr;


/**
 * \ingroup Exceptions
 *
 * \brief BadPerformanceCreation is thrown when \ref factory tries to generate an unknown performance evaluator
 */
class BadPerformanceCreation : public BadTaskCreation
{
public:

    /**
     * Exception constructor.
     */
    BadPerformanceCreation(std::string type): BadTaskCreation(type) {}
};

/**
 * \ingroup Performance
 * \brief Performance is the class that evaluates prediction performance
 */
template <typename T>
class Performance : public Task<T>
{
public:

	///
	/// \brief Constructor
	/// \param taskName The task name
	///
	Performance(const std::string& taskName)
    	:Task<T>("perf", taskName){}
	
	///
	/// \brief Default constructor
	///
	Performance()
    	:Task<T>("perf", ""){}

	/**
     * Evaluates prediction performance
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the different required fields based on the sub-class
     *
     * \return a GurlsOptionList
     */

    virtual GurlsOptionsList *execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) = 0;

    /**
     * Factory function returning a pointer to the newly created object.
     *
     * \warning The returned pointer is a plain, un-managed pointer. The calling
     * function is responsible of deallocating the object.
     */
    static Performance<T>* factory(const std::string& id) throw(BadPerformanceCreation)
    {
        if(id == "precrec")
            return new PerfPrecRec<T>;
        if(id == "macroavg")
            return new PerfMacroAvg<T>;
        if(id == "totalavg")
            return new PerfTotalAvg<T>;
        if(id == "rmse")
            return new PerfRmse<T>;
        if(id == "abserr")
            return new PerfAbsErr<T>;

        throw BadPerformanceCreation(id);
    }

};

}

#endif // _GURLS_PERF_H
