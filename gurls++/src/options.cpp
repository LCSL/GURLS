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
GURLS_EXPORT std::ostream& operator<<(std::ostream& os, GurlsOption& opt)
{
    opt.operator <<(os);
    return os;
}

GURLS_EXPORT std::ostream& OptString::operator<<(std::ostream& os)
{
    return os << this->getValue();
}

GURLS_EXPORT std::ostream& OptNumber::operator<<(std::ostream& os)
{
    return os << this->getValue();
}

GURLS_EXPORT std::ostream& OptStringList::operator<<(std::ostream& os)
{
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

GURLS_EXPORT std::ostream& OptNumberList::operator<<(std::ostream& os)
{
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

GURLS_EXPORT std::ostream& OptTaskSequence::operator<<(std::ostream& os)
{
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



GurlsOption::GurlsOption(OptTypes t):type(t) {}

OptTypes GurlsOption::getType() const
{
    return type;
}

GurlsOption::~GurlsOption(){}

bool GurlsOption::isA(OptTypes id) const
{
    return (id == GenericOption);
}

const std::type_info& GurlsOption::getDataID()
{
    return typeid(GurlsOption);
}

std::ostream& GurlsOption::operator<<(std::ostream& os)
{
    return os;
}



bool OptTaskSequence::isValid(const std::string & str, std::string& type, std::string& name)
{
    size_t found = str.find(gurls::TASKDESC_SEPARATOR);

    if (found==std::string::npos)
        return false;

    type = str.substr(0, found);
    name = str.substr(found+1);

    if (name.find(gurls::TASKDESC_SEPARATOR)!=std::string::npos)
        return false;

    return true;
}

OptTaskSequence::OptTaskSequence():GurlsOption(TaskSequenceOption)
{
    tasks = new std::vector<std::string>();
}

OptTaskSequence::OptTaskSequence(const char *str): GurlsOption(TaskSequenceOption)
{
    tasks = new std::vector<std::string>();
    tasks->push_back(str);
}

OptTaskSequence::OptTaskSequence(const std::string& str): GurlsOption(TaskSequenceOption)
{
    tasks = new std::vector<std::string>();
    tasks->push_back(str);
}

OptTaskSequence::OptTaskSequence(const std::vector<std::string> &data): GurlsOption(TaskSequenceOption)
{
    tasks = new std::vector<std::string>(data);
}

OptTaskSequence::~OptTaskSequence()
{
    tasks->clear();
    delete tasks;
}

void OptTaskSequence::addTask(const std::string newtask)
{
    tasks->push_back(newtask);
}

const std::vector<std::string>& OptTaskSequence::getValue() const
{
    return *tasks;
}

void OptTaskSequence::setValue(const std::vector<std::string> newvalue)
{
    tasks->clear();
    delete tasks;

    tasks = new std::vector<std::string>(newvalue.begin(), newvalue.end());
}

bool OptTaskSequence::isA(OptTypes id) const
{
    return (id == TaskSequenceOption);
}

OptTaskSequence *OptTaskSequence::dynacast(GurlsOption *opt)
{
    if (opt->isA(TaskSequenceOption) )
        return static_cast<OptTaskSequence*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

const OptTaskSequence *OptTaskSequence::dynacast(const GurlsOption *opt)
{
    if (opt->isA(TaskSequenceOption) )
        return static_cast<const OptTaskSequence*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

void OptTaskSequence::getTaskAt(int index, string &taskdesc, string &taskname)
{
    if (!isValid((*tasks)[index], taskdesc, taskname))
        throw new gException(gurls::Exception_Invalid_TaskSequence);
}

long OptTaskSequence::size()
{
    return tasks->size();
}





OptNumberList::OptNumberList(): GurlsOption(NumberListOption)
{
    value = new std::vector<double>();
}

OptNumberList::OptNumberList(const std::vector<double>& vec): GurlsOption(NumberListOption)
{
    value = new std::vector<double>(vec.begin(), vec.end());
}

OptNumberList::OptNumberList(double v): GurlsOption(NumberListOption)
{
    value = new std::vector<double>();
    value->push_back(v);
}

OptNumberList::OptNumberList(double *v, int n):GurlsOption(NumberListOption), value()
{
    value = new std::vector<double>(v, v+n);
}

OptNumberList::~OptNumberList()
{
    value->clear();

    delete value;
}

void OptNumberList::setValue(const std::vector<double> newvalue)
{
    value->clear();
    delete value;

    value = new std::vector<double>(newvalue.begin(), newvalue.end());
}

void OptNumberList::add(const double d)
{
    value->push_back(d);
}

std::vector<double> &OptNumberList::getValue()
{
    return *value;
}

const std::vector<double>& OptNumberList::getValue() const
{
    return *value;
}

bool OptNumberList::isA(OptTypes id) const
{
    return (id == NumberListOption);
}

OptNumberList *OptNumberList::dynacast(GurlsOption *opt)
{
    if (opt->isA(NumberListOption) )
       return static_cast<OptNumberList*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

const OptNumberList *OptNumberList::dynacast(const GurlsOption *opt)
{
    if (opt->isA(NumberListOption) )
        return static_cast<const OptNumberList*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}



OptNumber::OptNumber(): GurlsOption(NumberOption), value(0.0) {}

OptNumber::OptNumber(double v): GurlsOption(NumberOption), value(v){}

OptNumber& OptNumber::operator=(double other)
{
    this->type = NumberOption;
    this->value = other;
    return *this;
}

void OptNumber::setValue(double newvalue)
{
    value = newvalue;
}

const double& OptNumber::getValue() const
{
    return value;
}

double& OptNumber::getValue()
{
    return value;
}

bool OptNumber::isA(OptTypes id) const
{
    return (id == NumberOption);
}

OptNumber* OptNumber::dynacast(GurlsOption* opt)
{
    if (opt->isA(NumberOption) )
        return static_cast<OptNumber*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

const OptNumber* OptNumber::dynacast(const GurlsOption* opt)
{
    if (opt->isA(NumberOption) )
        return static_cast<const OptNumber*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}



OptStringList::OptStringList(): GurlsOption(StringListOption)
{
    value = new std::vector<std::string>();
}

OptStringList::OptStringList(const std::vector<std::string>& vec): GurlsOption(StringListOption)
{
    value = new std::vector<std::string>(vec.begin(), vec.end());
}

OptStringList::OptStringList(std::string& str): GurlsOption(StringListOption)
{
    value = new std::vector<std::string>();
    value->push_back(str);
}

OptStringList::~OptStringList()
{
    value->clear();
    delete value;
}

void OptStringList::setValue(const std::vector<std::string> newvalue)
{
    value->clear();
    delete value;

    value = new std::vector<std::string>(newvalue.begin(), newvalue.end());
}

void OptStringList::add(const std::string str)
{
    value->push_back(str);
}

const std::vector<std::string>& OptStringList::getValue() const
{
    return *value;
}

bool OptStringList::isA(OptTypes id) const
{
    return (id == StringListOption);
}

OptStringList* OptStringList::dynacast(GurlsOption* opt)
{
    if (opt->isA(StringListOption) )
        return static_cast<OptStringList*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}


const OptStringList* OptStringList::dynacast(const GurlsOption* opt)
{
    if (opt->isA(StringListOption) )
        return static_cast<const OptStringList*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}



OptString::OptString(): GurlsOption(StringOption), value(""){}

OptString::OptString(const char* str): GurlsOption(StringOption),value(str){}

OptString::OptString(const std::string& str): GurlsOption(StringOption),value(str){}

OptString::~OptString()
{
    value.clear();
}

OptString& OptString::operator=(const std::string& other)
{
    this->type = StringOption;
    this->value = other;
    return *this;
}

void OptString::setValue(const std::string& newvalue)
{
    value = newvalue;
}

std::string& OptString::getValue()
{
    return value;
}

const std::string& OptString::getValue() const
{
    return value;
}

bool OptString::isA(OptTypes id) const
{
    return (id == StringOption);
}

OptString* OptString::dynacast(GurlsOption* opt)
{
    if (opt->isA(StringOption) )
        return static_cast<OptString*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

const OptString* OptString::dynacast(const GurlsOption* opt)
{
    if (opt->isA(StringOption) )
        return static_cast<const OptString*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

}

