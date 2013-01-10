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

#ifndef _GURLS_OPTLIST_H_
#define _GURLS_OPTLIST_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <typeinfo>

#include "exports.h"
#include "gmat2d.h"
#include "gvec.h"
#include "options.h"
#include "optarray.h"
#include "exceptions.h"


namespace gurls {


#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4251)
#endif

/**
  * \ingroup Settings
  * \brief GurlsOptionsList is an option containing a list of options
  * mapped by name.
  */
class GURLS_EXPORT GurlsOptionsList: public GurlsOption
{
public:
    typedef std::map<std::string, GurlsOption* > ValueType;

    /**
      * Constructor. Builds an optionlist with a name and optionally a set of default options
      *
      * \param ExpName name of the options list
      * \param usedefopt if \a true the list is filled with a set of default options, if \a false the list is left empty
      */
    GurlsOptionsList(std::string ExpName, bool usedefopt = false);

    /**
      * Destructor
      */
    ~GurlsOptionsList();

    /**
      * Adds a generic option to the list indexed with a specified key
      */
    bool addOpt(std::string key, GurlsOption* value);

    /**
      * Adds a string option to the list indexed with a specified key
      */
    bool addOpt(std::string key, std::string value);

    /**
      * Returns a pointer to a generic option mapped with a key
      */
    GurlsOption* getOpt(std::string key);

    /**
      * Returns a pointer to a generic option mapped with a key
      */
    const GurlsOption* getOpt(std::string key) const;

    /**
      * Returns a pointer to a T option mapped with a key
      */
    template<class T>
    T* getOptAs(std::string key)
    {
        return T::dynacast(this->getOpt(key));
    }

    /**
      * Returns a pointer to a T option mapped with a key
      */
    template<class T>
    const T* getOptAs(std::string key) const
    {
        return T::dynacast(this->getOpt(key));
    }

    /**
      * Returns a reference to the value contained into an option mapped with a key
      */
    template<class T>
    typename T::ValueType& getOptValue(std::string key)
    {
        return this->getOptAs<T>(key)->getValue();
    }

    /**
      * Returns a reference to the value contained into an option mapped with a key
      */
    template<class T>
    const typename T::ValueType& getOptValue(std::string key) const
    {
        return this->getOptAs<T>(key)->getValue();
    }

    /**
      * Returns a string option mapped with a key
      */
    std::string getOptAsString(std::string key) const;

    /**
      * Returns the list name
      */
    std::string getName() const;

    /**
      * Returns the entire map
      */
    const ValueType& getValue() const;

    /**
      * Sets the list name
      */
    void setName(std::string);

    /**
      * Returns a numeric option mapped with a key
      */
    double getOptAsNumber(std::string key) const;

    /**
      * Prints the options list
      */
    void printAll();

    /**
      * Checks if the list has an option mapped with a specified key
      */
    bool hasOpt(std::string key) const;

    /**
      * Removes the option mapped with a specified key
      *
      * \param key string key
      * \param deleteMembers If true deallocates the removed option, if false option will be only detached
      */
    void removeOpt(std::string key, bool deleteMembers = true);

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref GurlsOptionsList
      */
    static GurlsOptionsList* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref GurlsOptionsList
      */
    static const GurlsOptionsList* dynacast(const GurlsOption* opt);

    /**
      * Returns the number of options within the list
      */
    int size() const;

    /**
      * Returns a pointer to the idx-th option into the list
      */
    GurlsOption* operator[] (int idx);

    /**
      * Writes a GurlsOptionsList to a stream
      */
    friend GURLS_EXPORT std::ostream& operator<<(std::ostream& os, GurlsOptionsList& opt);

    /**
      * Writes the list to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    /**
      * Writes the list to a string
      */
    std::string toString();

    /**
      *
      */
    void copyOpt(std::string key, const GurlsOptionsList &from);

    /**
      * Serializes the list to file
      */
    void save(const std::string& fileName) const;

    /**
      * Deserializes the list from file
      */
    void load(const std::string& fileName);

protected:
    std::string name;   ///< Option name
    ValueType* table;   ///< Options list, indexed by name

};

#ifdef _WIN32
#pragma warning(pop)
#endif

}

#endif // _GURLS_OPTLIST_H_
