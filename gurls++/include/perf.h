/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011, IIT@MIT Lab
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
#include <exception>
#include <stdexcept>

#include "gmath.h"
#include "options.h"
#include "optlist.h"

using namespace std;

namespace gurls {

template <typename T>
class MacroAvg;

template <typename T>
class PrecisionRecall;


//template <typename Matrix>
template <typename T>
class Performance
{
public:
//	virtual Matrix& execute( const Matrix& X, const Matrix& Y, GurlsOptionsList& opt) = 0;
    virtual void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt) = 0;

    class BadPerformanceCreation : public std::logic_error {
    public:
        BadPerformanceCreation(std::string type)
            : logic_error("Cannot create type " + type) {}
    };
//    static Performance<Matrix>*
    static Performance<T>*
    factory(const std::string& id) throw(BadPerformanceCreation)
    {
        if(id == "precrec")
            return new PrecisionRecall<T>;
        else if(id == "macroavg")
            return new MacroAvg<T>;
        else
            throw BadPerformanceCreation(id);
    }

};

}

#endif // _GURLS_PERF_H
