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

#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>

#include "options.h"

using namespace std;

namespace gurls{

/**
  * Writes a GurlsOption to a stream
  */
GURLS_EXPORT std::ostream& operator<<(std::ostream& os, GurlsOption& opt) {
    opt.operator <<(os);
    return os;
}



/*
std::ostream& operator<<(std::ostream& os, OptString& opt) {
    return os << opt.getValue();
}

std::ostream& operator<<(std::ostream& os, OptNumber& opt) {
    return os << static_cast<double>(opt.getValue());
}

std::ostream& operator<<(std::ostream& os, OptStringList& opt){
    std::vector<std::string> V = opt.getValue();
    std::vector<std::string>::iterator it = V.begin();
    std::vector<std::string>::iterator end = V.end();
    os << "(" << (*it++);
    while( it != end){
        os << ", " << (*it++);
    }
    os << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, OptNumberList& opt){
    std::vector<double> V = opt.getValue();
    std::vector<double>::iterator it = V.begin();
    std::vector<double>::iterator end = V.end();
    os << "(" << (*it++);
    while( it != end){
        os << ", " << (*it++);
    }
    os << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, OptFunction& opt){
    os << "Pointer to the function " << opt.getName() << " whose signature is: double (*func)(double*, int)" ;
    return os;
}
*/


GURLS_EXPORT std::ostream& OptString::operator<<(std::ostream& os){
    return os << this->getValue();
}

GURLS_EXPORT std::ostream& OptNumber::operator<<(std::ostream& os) {
    return os << this->getValue();
}

GURLS_EXPORT std::ostream& OptStringList::operator<<(std::ostream& os){
    const std::vector<std::string>& V = this->getValue();
    std::vector<std::string>::const_iterator it = V.begin();
    std::vector<std::string>::const_iterator end = V.end();

    os << "(";
    if(!V.empty())
        os << (*it++);

    while( it != end)
        os << ", " << (*it++);
    os << ")";

    return os;
}

GURLS_EXPORT std::ostream& OptNumberList::operator<<(std::ostream& os){
    std::vector<double>& V = this->getValue();
    std::vector<double>::iterator it = V.begin();
    std::vector<double>::iterator end = V.end();

    os << "(";
    if(!V.empty())
        os << (*it++);

    while( it != end)
        os << ", " << (*it++);
    os << ")";

    return os;
}

GURLS_EXPORT std::ostream& OptFunction::operator<<(std::ostream& os){
    os << "Pointer to the function <" << this->getName()
       << "> whose signature is: T (*func)(T*, int)" ;
    return os;
}

GURLS_EXPORT std::ostream& OptTaskSequence::operator<<(std::ostream& os){

    std::vector<std::string>::iterator it = tasks->begin();
    std::vector<std::string>::iterator end = tasks->end();

    os << "(";
    if(!tasks->empty())
        os << (*it++);

    while( it != end)
        os << ", " << (*it++);
    os << ")";

    return os;
}

bool OptTaskSequence::isValid(const std::string & str, std::string& type, std::string& name) {
    size_t found = str.find(gurls::TASKDESC_SEPARATOR);
    if (found==std::string::npos){
        return false;
    }
    type = str.substr(0, found);
    name = str.substr(found+1);
    if (name.find(gurls::TASKDESC_SEPARATOR)!=std::string::npos){
        return false;
    }
    return true;
}

/**
  * OptMatrix empty constructor for float elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<float> >::OptMatrix(): OptMatrixBase () , value (*(new gMat2D<float>(2,2)))
{
    this->matType = FLOAT;
}

/**
  * OptMatrix constructor for float elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<float> >::OptMatrix(gMat2D<float>& m): OptMatrixBase(), value(m)
{
    this->matType = FLOAT;
}

/**
  * OptMatrix empty constructor for double elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<double> >::OptMatrix(): OptMatrixBase () , value (*(new gMat2D<double>(2,2)))
{
    this->matType = DOUBLE;
}

/**
  * OptMatrix constructor for double elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<double> >::OptMatrix(gMat2D<double>& m): OptMatrixBase(), value(m)
{
    this->matType = DOUBLE;
}

/**
  * OptMatrix empty constructor for unsigned long elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<unsigned long> >::OptMatrix(): OptMatrixBase () , value (*(new gMat2D<unsigned long>(2,2)))
{
    this->matType = ULONG;
}

/**
  * OptMatrix constructor for unsigned long elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<unsigned long> >::OptMatrix(gMat2D<unsigned long>& m): OptMatrixBase(), value(m)
{
    this->matType = ULONG;
}

}

