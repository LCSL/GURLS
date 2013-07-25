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

#ifndef _GURLS_OPTARRAY_H_
#define _GURLS_OPTARRAY_H_

#include "gurls++/options.h"

#include <vector>

namespace gurls
{

class GurlsOptionsList;

/**
  * \ingroup Settings
  * \brief Optarray is an option containing an indexed array of options
  */
class GURLS_EXPORT OptArray: public GurlsOption
{
public:

    typedef std::vector<GurlsOption*> ValueType; ///<

protected:

    ValueType *value; ///< Options array


public:

    /**
      * Empty constructor, builds an empty array.
      */
    OptArray();

    /**
      * Constructor form an existing \ref OptArray
      */
    OptArray(const OptArray& other);

    /**
      * Destructor
      */
    ~OptArray();

    /**
      * Returns the array
      */
    const ValueType& getValue() const;

    /**
      * Adds a new option at the end of the array
      */
    void push_back(GurlsOption *opt);

    /**
      * Removes the option at a specified index
      *
      * \param index Index
      * \param deleteMembers If true deallocates the removed option, if false option will be only detached
      */
    void erase(unsigned long index, bool deleteMembers = true);

    /**
      * Removes all elements from the array
      */
    void clear();

    /**
      * Request a change in capacity
      */
    void reserve(unsigned long size);

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptArray
      */
    static OptArray* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptArray
      */
    static const OptArray* dynacast(const GurlsOption* opt);

    /**
      * Returns the number of options within the array
      */
    unsigned long size() const;

    /**
      * Returns a pointer to the i-th option into the array
      */
    GurlsOption* operator[] (unsigned long i) const;

    /**
      * Writes an OptArray to a stream
      */
    friend GURLS_EXPORT std::ostream& operator<<(std::ostream& os, const OptArray& opt);

    /**
      * Writes the array to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os) const;

    /**
      * Serializes the array to file
      */
    void save(const std::string& fileName) const;

    /**
      * Deserializes the array from file
      */
    void load(const std::string& fileName);

};

}

#endif //_GURLS_OPTARRAY_H_
