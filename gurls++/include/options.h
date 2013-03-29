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

#ifndef _GURLS_OPT_H_
#define _GURLS_OPT_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <typeinfo>

#include "exports.h"
#include "gmat2d.h"
#include "gvec.h"
#include "exceptions.h"


/**
  * \ingroup Settings
  * \file
  */
namespace gurls
{


/**
  * \enum OptTypes
  * Enumeration containing all implemented Option types
  */
enum OptTypes	{GenericOption, StringOption, NumberOption,
                 StringListOption, NumberListOption, FunctionOption,
                 MatrixOption, VectorOption,
                OptListOption, TaskSequenceOption, TaskIDOption, OptArrayOption, ProcessOption};



/**
  * \ingroup Settings
  * \brief GurlsOption is an abstraction of a generic `option', which is
  * widely used within the GURLS++ package to store either numeric
  * parameters necessary to configure specific algorigms or sequences
  * of strings holding the names of the specific procedures that
  * have to be performed.
  *
  * The instances of the GurlsOption class hold information about
  * the type (one of the elements in the OptTypes enumeration),
  * while the value related to each specific option is stored using
  * the attributes of the subclasses of GurlsOption.
  */
class GURLS_EXPORT GurlsOption
{
protected:
    OptTypes type; ///< Option type

public:

    /**
      * Constructor from an option type
      */
    GurlsOption(OptTypes t);

    /**
      * Returns the option type
      */
    OptTypes getType() const;

    /**
      * Destructor
      */
    virtual ~GurlsOption();

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Returns the identifier of the option type
      */
    const virtual std::type_info& getDataID();

    /**
      * Writes an option to a stream
      */
    friend GURLS_EXPORT std::ostream& operator<<(std::ostream& os, GurlsOption& opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);
};

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4251)
#endif

/**
  * \ingroup Settings
  * \brief OptString is an option containing a generic string.
  */
class GURLS_EXPORT OptString: public GurlsOption
{

protected:
    std::string value; ///< String containing the option value

public:

    typedef std::string ValueType;

    /**
      * Empty constructor
      */
    OptString();

    /**
      * Constructor from a char buffer
      */
//    OptString(const char* str);

    /**
      * Constructor from a string
      */
    OptString(const std::string& str);

	/**
      * Constructor from a wstring
      */
	OptString(const std::wstring& str);

    /**
      * Copies the opt values from an existing \ref OptString
      */
    OptString& operator=(const OptString& other);

    /**
      * Destructor
      */
    ~OptString();

    /**
      * Copies the string value of the option from an existing string
      */
    OptString& operator=(const std::string& other);

    /**
      * Copies the string value of the option from an existing string
      */
    void setValue(const std::string& newvalue);

    /**
      * Returns the string value
      */
    std::string& getValue();

    /**
      * Returns the string value
      */
    const std::string& getValue() const;

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptString
      */
    static OptString* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptString
      */
    static const OptString* dynacast(const GurlsOption* opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

};

#ifdef _WIN32
#pragma warning(pop)
#endif

/**
  * \ingroup Settings
  * \brief OptStringList is an option containing a list of strings.
  */
class GURLS_EXPORT OptStringList: public GurlsOption
{
protected:
    std::vector<std::string>* value; ///< Vector of strings containing the option value

public:

    typedef std::vector<std::string> ValueType;

    /**
      * Empty constructor
      */
    OptStringList();

    /**
      * Constructor from a vector of strings
      */
    OptStringList(const std::vector<std::string>& vec);

    /**
      * Constructor from a string, builds a 1-size vector of strings
      */
    OptStringList(std::string& str);

    /**
      * Copies the option values from an existing \ref OptStringList
      */
    OptStringList& operator=(const OptStringList& other);

    /**
      * Destructor
      */
    ~OptStringList();

    /**
      * Copies the opt values from a vector of strings
      */
    void setValue(const std::vector<std::string> newvalue);

    /**
      * Adds a string to the list
      */
    void add(const std::string str);

    /**
      * Returns the vector of strings
      */
    const std::vector<std::string>& getValue() const;

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptStringList
      */
    static OptStringList* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptStringList
      */
    static const OptStringList* dynacast(const GurlsOption* opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    /**
      * Removes all contents in the list
      */
    void clear();

    /**
      * Adds a string to the list
      */
    OptStringList& operator<<(std::string& str);

    /**
      * Adds a string to the list
      */
    OptStringList& operator<<(const char* str);

};


/**
  * \ingroup Settings
  * \brief OptNumber is an option containing a double precision
  * floating point number
  */
class GURLS_EXPORT OptNumber: public GurlsOption
{
protected:
    double value; ///< Option value

public:

    typedef double ValueType;

    /**
      * Empty constructor
      */
    OptNumber();

    /**
      * Constructor from a double
      */
    OptNumber(double v);

    /**
      * Copies the option values from an existing \ref OptNumber
      */
    OptNumber& operator=(const OptNumber& other);

//    ~OptNumber(){}

    /**
      * Copies the option value from a double
      */
    OptNumber& operator=(double other);

    /**
      * Sets the option value to the given one
      */
    void setValue(double newvalue);

    /**
      * Returns the option value
      */
    const double& getValue() const;

    /**
      * Returns the option value
      */
    double& getValue();

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumber
      */
    static OptNumber* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumber
      */
    static const OptNumber* dynacast(const GurlsOption* opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

};

//#ifdef GURLS_DEPRECATED

/**
  * \ingroup Settings
  * \brief OptNumberList is an option containing a list of
  * double precision floating point numbers
  */
class GURLS_EXPORT OptNumberList: public GurlsOption
{
private:
    std::vector<double>* value; ///< Option value

public:

    typedef std::vector<double> ValueType;

    /**
      * Empty constructor
      */
    OptNumberList();

    /**
      * Constructor from a vector of double
      */
    OptNumberList(const std::vector<double>& vec);

    /**
      * Constructor from a double, builds a 1-size vector of double
      */
    OptNumberList(double v);

    /**
      * Constructor from a double buffer of size n
      */
    OptNumberList(double *v, int n);

    /**
      * Copies the option values from an existing \ref OptNumberList
      */
    OptNumberList& operator=(const OptNumberList& other);

    ~OptNumberList();

    /**
      * Copies the option values from a vector of double
      */
    void setValue(const std::vector<double> newvalue);

    /**
      * Adds a double to the list
      */
    void add(const double d);

    /**
      * Returns the vector of double inside the option
      */
    std::vector<double>& getValue();

    /**
      * Returns the vector of double inside the option
      */
    const std::vector<double>& getValue() const;

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumberList
      */
    static OptNumberList* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumberList
      */
    static const OptNumberList* dynacast(const GurlsOption* opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    /**
      * Removes all contents in the list
      */
    void clear();

    /**
      * Adds a string to the list
      */
    OptNumberList& operator<<(double& d);
}

//#ifdef __GNUC__
//    __attribute__ ((deprecated))
//#endif
;

//#ifdef _WIN32
//#pragma deprecated( OptNumberList )
//#endif

//#endif


/**
  * String used to tokenize task strings (e.g. "<task_desc>TASKDESC_SEPARATOR<task_name>")
  */
static const std::string TASKDESC_SEPARATOR(":");

/**
  * \ingroup Settings
  * \brief OptTaskSequence is an option containing
  * a sequence of task that forms a pipeline
  */
class GURLS_EXPORT OptTaskSequence: public OptStringList
{
protected:
    /**
      * Parses a string cheching if it's a valid task string, in the form "<task_desc>TASKDESC_SEPARATOR<task_name>"
      */
    bool isValid(const std::string & str, std::string& type, std::string& name);

public:

    typedef OptStringList::ValueType ValueType;

    /**
      * Empty constructor
      */
    OptTaskSequence();

    /**
      * Constructor from a buffer of chars, builds a 1-size vector of strings
      */
    OptTaskSequence(const char* str);

    /**
      * Constructor from a string, builds a 1-size vector of strings
      */
    OptTaskSequence(std::string &str);

    /**
      * Constructor from a string vector
      */
    OptTaskSequence(const std::vector<std::string>& data);

    /**
      * Copies the matrix from an existing \ref OptTaskSequence
      */
    OptTaskSequence& operator=(const OptTaskSequence& other);

    /**
      * Adds a new task string to the sequence
      */
    void addTask(const std::string newtask);

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptTaskSequence
      */
    static OptTaskSequence* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptTaskSequence
      */
    static const OptTaskSequence* dynacast(const GurlsOption* opt);

    /**
      * Parse the task string at a given index and returns the task description and the task name
      */
    void getTaskAt(int index, std::string& taskdesc, std::string& taskname);

    /**
      * Returns the number of tasks into the sequence
      */
    unsigned long size();

};

/**
  * \ingroup Settings
  * \brief OptProcess is an option containing
  * a sequence of actions that form a Gurls process
  */
class GURLS_EXPORT OptProcess: public GurlsOption
{
public:
    /**
     * Execution options for a GURLS task
     */
    enum Action {ignore, compute, computeNsave, load, remove};

    typedef std::vector<Action> ValueType;

protected:

    ValueType* value;   ///< Option value

    /**
     * Used to map enum values to their names
     */
    static std::vector<std::string>& actionNames();

public:

    /**
      * Empty constructor
      */
    OptProcess();

    /**
      * Constructor from an existing OptProcess
      */
    OptProcess(const OptProcess& other);

    /**
      * Destructor
      */
    ~OptProcess();

    /**
      * Adds a new action into the process
      */
    void addAction(const Action action);

    /**
      * Adds a new action into the process
      */
    OptProcess& operator<<(const Action action);

    /**
      * Returns the option's value
      */
    const ValueType& getValue() const;

    /**
      * Gets the element at a given position in the sequence
      */
    Action operator[](unsigned long index);

    /**
      * Removes all actions
      */
    void clear();

    /**
      * Returns the number of actions into the process
      */
    unsigned long size();

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptProcess
      */
    static OptProcess* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptProcess
      */
    static const OptProcess* dynacast(const GurlsOption* opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

};

}

#endif // _GURLS_OPT_H_
