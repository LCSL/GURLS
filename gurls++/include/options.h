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

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>


/**
  * \ingroup Settings
  * \file
  */

namespace gurls {


/**
  * \enum OptTypes
  * Enumeration containing all implemented Option types
  */
enum OptTypes	{GenericOption, StringOption, NumberOption,
                 StringListOption, NumberListOption, FunctionOption,
                 MatrixOption, VectorOption,
                OptListOption, TaskSequenceOption, TaskIDOption};


/*
 * Typical functions used within the gurls package to combine or choose
 * among different regularization parameters or multiple numeric
 * options.
 */

/**
 * Computes the mean value of a vector \c v of lenght \c n
 */
template <typename T>
T mean(const T *v, int n)
{
    T m = (T)0.0;

    for (int i = 0; i < n; ++i, ++v)
        m += *v;

    return m/n;
}

/**
 * Computes the smallest element in a vector \c v of lenght \c n
 */
template <typename T>
T min(const T *v, int n)
{
    return (*std::min_element(v,v+n));
}

/**
 * Computes the largest element in a vector \c v of lenght \c n
 */
template <typename T>
T max(const T *v, int n)
{
    return (*std::max_element(v,v+n));
}

/**
 * Computes the median value of a vector \c v of lenght \c n
 */
template <typename T>
T median(const T *v, int n)
{
    std::vector<T> vd(v, v+n);
    sort(vd.begin(), vd.end());

    if(n%2)
        return *(vd.begin()+vd.size()/2);
    else
        return (*(vd.begin()+vd.size()/2) + *((vd.begin()+vd.size()/2)-1) )/2;
}


/**
  * \ingroup Settings
  * \brief GurlsOption is an abstraction of a generic `option', which is
  * widely used within the GURLS++ package to store either numeric
  * parameters necessary to configure specific algorigms or sequences
  * of strings holding the names of the specific procedures that
  * have to be performed.
  *
  * The instances of the GURLSOPTION class hold information about
  * the type (one of the elements in the OptTypes enumeration),
  * while the value related to each specific option is stored using
  * the attributes of the subclasses of GURLSOPTION.
  */
class GURLS_EXPORT GurlsOption
{
protected:
    OptTypes type; ///< Option type

public:

    /**
      * Constructor from an option type
      */
    GurlsOption(OptTypes t):type(t) {}

    /**
      * Returns the option type
      */
    const OptTypes getType() const {return type;}

    /**
      * Destructor
      */
    virtual ~GurlsOption(){}

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == GenericOption); }

    /**
      * Returns the identifier of the option type
      */
    const virtual std::type_info& getDataID(){
        return typeid(GurlsOption);
    }

    /**
      * Writes an option to a stream
      */
    friend GURLS_EXPORT std::ostream& operator<<(std::ostream& os, GurlsOption& opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os){return os;}
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
private:
    std::string value; ///< String containing the option value

public:

    typedef std::string ValueType;

    /**
      * Empty constructor
      */
    OptString(): GurlsOption(StringOption), value(""){}

    /**
      * Constructor from a char buffer
      */
    OptString(const char* str): GurlsOption(StringOption),value(str){}

    /**
      * Constructor from a string
      */
    OptString(const std::string& str): GurlsOption(StringOption),value(str){}

    /**
      * Copies the opt values from an existing \ref OptString
      */
    OptString& operator=(const OptString& other);

    /**
      * Destructor
      */
    ~OptString(){value.clear();}

    /**
      * Copies the string value of the option from an existing string
      */
    OptString& operator=(const std::string& other){
        this->type = StringOption;
        this->value = other;
        return *this;
    }

    /**
      * Copies the string value of the option from an existing string
      */
    void setValue(const std::string& newvalue) {value = newvalue;}

    /**
      * Returns the string value
      */
    std::string& getValue() { return value;}

    /**
      * Returns the string value
      */
    const std::string& getValue() const { return value;}

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == StringOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptString
      */
    static OptString* dynacast(GurlsOption* opt) {
        if (opt->isA(StringOption) ){
            return static_cast<OptString*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptString
      */
    static const OptString* dynacast(const GurlsOption* opt) {
        if (opt->isA(StringOption) ){
            return static_cast<const OptString*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        ar & this->value;
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        ar & this->value;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
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
private:
    std::vector<std::string>* value; ///< Vector of strings containing the option value

public:

    typedef std::vector<std::string> ValueType;

    /**
      * Empty constructor
      */
    OptStringList(): GurlsOption(StringListOption){
        value = new std::vector<std::string>();
    }

    /**
      * Constructor from a vector of strings
      */
    OptStringList(const std::vector<std::string>& vec): GurlsOption(StringListOption){
        value = new std::vector<std::string>(vec.begin(), vec.end());
    }

    /**
      * Constructor from a string, builds a 1-size vector of strings
      */
    OptStringList(std::string& str): GurlsOption(StringListOption){
        value = new std::vector<std::string>();
        value->push_back(str);
    }

    /**
      * Copies the option values from an existing \ref OptStringList
      */
    OptStringList& operator=(const OptStringList& other);

    /**
      * Destructor
      */
    ~OptStringList(){
        value->clear();
        delete value;
    }

    /**
      * Copies the opt values from a vector of strings
      */
    void setValue(const std::vector<std::string> newvalue) {
        delete value;
        value = new std::vector<std::string>(newvalue.begin(), newvalue.end());
    }

    /**
      * Adds a string to the list
      */
    void add(const std::string str){
        value->push_back(str);
    }

    /**
      * Returns the vector of strings
      */
    const std::vector<std::string>& getValue() const { return *value;}

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == StringListOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptStringList
      */
    static OptStringList* dynacast(GurlsOption* opt) {
        if (opt->isA(StringListOption) ){
            return static_cast<OptStringList*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptStringList
      */
    static const OptStringList* dynacast(const GurlsOption* opt) {
        if (opt->isA(StringListOption) ){
            return static_cast<const OptStringList*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        int n = this->value->size();
        ar & n;
        for (int i = 0; i < n; i++){
            ar & this->value->at(i);
        }
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        int n = 0;
        ar & n;
        value->clear();
        for(int i = 0; i < n; i++){
            std::string s("");
            ar & s;
            this->value->push_back(s);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

};


/**
  * \ingroup Settings
  * \brief OptNumber is an option containing a double precision
  * floating point number
  */
class GURLS_EXPORT OptNumber: public GurlsOption
{
private:
    double value; ///< Option value

public:

    typedef double ValueType;

    /**
      * Empty constructor
      */
    OptNumber(): GurlsOption(NumberOption), value(0) {}

    /**
      * Constructor from a double
      */
    OptNumber(double v): GurlsOption(NumberOption), value(v) {}

    /**
      * Copies the option values from an existing \ref OptNumber
      */
    OptNumber& operator=(const OptNumber& other);

//    ~OptNumber(){}

    /**
      * Copies the option value from a double
      */
    OptNumber& operator=(double other){
        this->type = NumberOption;
        this->value = other;
        return *this;
    }

    /**
      * Sets the option value to the given one
      */
    void setValue(double newvalue) {value = newvalue;}

    /**
      * Returns the option value
      */
    const double& getValue() const {return value;}

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == NumberOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumber
      */
    static OptNumber* dynacast(GurlsOption* opt) {
        if (opt->isA(NumberOption) ){
            return static_cast<OptNumber*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumber
      */
    static const OptNumber* dynacast(const GurlsOption* opt) {
        if (opt->isA(NumberOption) ){
            return static_cast<const OptNumber*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        ar & this->value;
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        ar & this->value;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

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
    OptNumberList(): GurlsOption(NumberListOption){
        value = new std::vector<double>();
    }

    /**
      * Constructor from a vector of double
      */
    OptNumberList(const std::vector<double>& vec): GurlsOption(NumberListOption){
        value = new std::vector<double>(vec.begin(), vec.end());
    }

    /**
      * Constructor from a double, builds a 1-size vector of double
      */
    OptNumberList(double v): GurlsOption(NumberListOption){
        value = new std::vector<double>();
        value->push_back(v);
    }

    /**
      * Constructor from a double buffer of size n
      */
    OptNumberList(double *v, int n):GurlsOption(NumberListOption), value(){
        value = new std::vector<double>(v, v+n);
    }

    /**
      * Copies the option values from an existing \ref OptNumberList
      */
    OptNumberList& operator=(const OptNumberList& other);

    ~OptNumberList()
    {
        delete value;
    }

    /**
      * Copies the option values from a vector of double
      */
    void setValue(const std::vector<double> newvalue) {
        value = new std::vector<double>(newvalue.begin(), newvalue.end());
    }

    /**
      * Adds a double to the list
      */
    void add(const double d){
        value->push_back(d);
    }

    /**
      * Returns the vector of double inside the option
      */
    std::vector<double>& getValue() { return *value;}

    /**
      * Returns the vector of double inside the option
      */
    const std::vector<double>& getValue() const { return *value;}

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == NumberListOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumberList
      */
    static OptNumberList* dynacast(GurlsOption* opt) {
        if (opt->isA(NumberListOption) ){
            return static_cast<OptNumberList*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptNumberList
      */
    static const OptNumberList* dynacast(const GurlsOption* opt)
    {
        if (opt->isA(NumberListOption) ){
            return static_cast<const OptNumberList*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        int n = value->size();
        ar & n;
        for (int i = 0; i < n; i++){
            ar & value->at(i);
        }
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        int n = 0;
        ar & n;
        value->clear();
        for(int i = 0; i < n; i++){
            double s;
            ar & s;
            value->push_back(s);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

}
//#ifdef __GNUC__
//    __attribute__ ((deprecated))
//#endif
;

//#ifdef _WIN32
//#pragma deprecated( OptNumberList )
//#endif

//#endif


#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4251)
#endif

/**
  * \ingroup Settings
  * \brief OptFunction is an option representing a pointer to a generic function
  * double (*function)(double* , int) operating over an array of floating point numbers.
  */
class GURLS_EXPORT OptFunction: public GurlsOption
{
private:
    std::string name; ///< Function name

public:
    /**
      * Constructor from a fuction name
      */
    OptFunction(std::string func_name): GurlsOption(FunctionOption), name(func_name) {}

    /**
      * Copies the option values from an existing \ref OptFunction
      */
    OptFunction& operator=(const OptFunction& other);

    /**
      * Copies the option values from a string representing a function name
      */
    void setValue(std::string func_name) {
        name = func_name;
    }

    /**
      * Returns the function name
      */
    std::string getName() const {return name;}

    /**
      * Executes the function over a buffer of length n, returning the result
      */
    template<typename T>
    T getValue(T* array, int n)
    {
        if(!name.compare("mean"))
            return mean(array, n);
        if(!name.compare("min"))
            return min(array, n);
        if(!name.compare("max"))
            return max(array, n);
        if(!name.compare("median"))
            return median(array, n);

        throw gException(gurls::Exception_Unknown_Function);
    }

    /**
      * Executes the function over a buffer of length n, returning the result
      */
    template<typename T>
    T getValue(T* array, int n) const
    {
        if(!name.compare("mean"))
            return mean(array, n);
        if(!name.compare("min"))
            return min(array, n);
        if(!name.compare("max"))
            return max(array, n);
        if(!name.compare("median"))
            return median(array, n);

        throw gException(gurls::Exception_Unknown_Function);
    }

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == FunctionOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptFunction
      */
    static OptFunction* dynacast(GurlsOption* opt) {
        if (opt->isA(FunctionOption) ){
            return static_cast<OptFunction*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptFunction
      */
    static const OptFunction* dynacast(const GurlsOption* opt)
    {
        if (opt->isA(FunctionOption) ){
            return static_cast<const OptFunction*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        ar & this->name;
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        ar & this->name;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
};


#ifdef _WIN32
#pragma warning(pop)
#endif

/**
  * \ingroup Settings
  * \brief OptMatrixBase is the base class for all options containing matrices.
  */
class GURLS_EXPORT OptMatrixBase: public GurlsOption
{
public:

    /**
      * Empty constructor
      */
    OptMatrixBase(): GurlsOption(MatrixOption){}

    /**
      * \enum MatrixType
      * Enumeration containing all supported element types
      */
    enum MatrixType{FLOAT, DOUBLE, ULONG};

    /**
      * Returns the element type for the matrix
      */
    MatrixType getMatrixType()
    {
        return matType;
    }

protected:
    MatrixType matType; ///< Stores the type of the elements inside the matrix
};

/**
  * \ingroup Settings
  * \brief OptMatrix is an option containing a matrix.
  * \tparam MatrixType Type of the matrix contained into the option
  */
template <typename Matrix>
class OptMatrix: public OptMatrixBase
{
private:
    Matrix& value;  ///< Option value

public:
    typedef Matrix ValueType;

    /**
      * Empty constructor
      */
    OptMatrix(): OptMatrixBase () , value (*(new Matrix(2,2)))
    {
        throw gException(Exception_Unsupported_MatrixType);
    }

    /**
      * Constructor from an existing matrix
      */
    OptMatrix(Matrix& m): OptMatrixBase(), value(m)
    {
        throw gException(Exception_Unsupported_MatrixType);
    }

    /**
      * Copies the option values from an existing \ref OptMatrix
      */
    OptMatrix<Matrix>& operator=(const OptMatrix<Matrix>& other);

    /**
      * Destructor
      */
    ~OptMatrix()
    {
        delete &value;
    }

    /**
      * Copies the matrix from an existing matrix
      */
    OptMatrix& operator=(const Matrix& other){
        this->type = MatrixOption;
        this->value = other;
        return *this;
    }

    /**
      * Copies the matrix from an existing matrix
      */
    void setValue(const Matrix& newvalue) {value = newvalue;}

    /**
      * Returns the matrix
      */
    Matrix& getValue() { return value;}

    /**
      * Returns the matrix
      */
    const Matrix& getValue() const { return value;}

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == MatrixOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptMatrix
      */
    static OptMatrix* dynacast(GurlsOption* opt) {
        if (opt->isA(MatrixOption) ){
            return static_cast<OptMatrix*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptMatrix
      */
    static const OptMatrix* dynacast(const GurlsOption* opt) {
        if (opt->isA(MatrixOption) ){
            return static_cast<const OptMatrix*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        ar & this->value;
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        ar & this->value;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

/**
  * OptMatrix empty constructor for float elements
  */
template <>
OptMatrix <gMat2D<float> >::OptMatrix();

/**
  * OptMatrix constructor for float elements
  */
template <>
OptMatrix <gMat2D<float> >::OptMatrix(gMat2D<float>& m);

/**
  * OptMatrix empty constructor for double elements
  */
template <>
OptMatrix <gMat2D<double> >::OptMatrix();

/**
  * OptMatrix constructor for double elements
  */
template <>
OptMatrix <gMat2D<double> >::OptMatrix(gMat2D<double>& m);

/**
  * OptMatrix empty constructor for unsigned long elements
  */
template <>
OptMatrix <gMat2D<unsigned long> >::OptMatrix();

/**
  * OptMatrix constructor for unsigned long elements
  */
template <>
OptMatrix <gMat2D<unsigned long> >::OptMatrix(gMat2D<unsigned long>& m);

/**
  * String used to tokenize task strings (e.g. "<task_desc>TASKDESC_SEPARATOR<task_name>")
  */
static const std::string TASKDESC_SEPARATOR(":");

/**
  * \ingroup Settings
  * \brief OptTaskSequence is an option containing
  * a sequence of task that forms a pipeline
  */
class GURLS_EXPORT OptTaskSequence: public GurlsOption
{
private:
    std::vector<std::string>* tasks; ///< Vector of tasks

    /**
      * Parses a string cheching if it's a valid task string, in the form "<task_desc>TASKDESC_SEPARATOR<task_name>"
      */
    bool isValid(const std::string & str, std::string& type, std::string& name);

public:

    typedef std::vector<std::string> ValueType;

    /**
      * Empty constructor
      */
    OptTaskSequence(): GurlsOption(TaskSequenceOption){
        tasks = new std::vector<std::string>();
    }

    /**
      * Constructor from a buffer of chars, builds a 1-size vector of strings
      */
    OptTaskSequence(const char* str): GurlsOption(TaskSequenceOption){
        tasks = new std::vector<std::string>();
        tasks->push_back(str);
    }

    /**
      * Constructor from a string, builds a 1-size vector of strings
      */
    OptTaskSequence(const std::string& str): GurlsOption(TaskSequenceOption){
        tasks = new std::vector<std::string>();
        tasks->push_back(str);
    }

    /**
      * Constructor from a string vector
      */
    OptTaskSequence(const std::vector<std::string>& data): GurlsOption(TaskSequenceOption)
    {
        tasks = new std::vector<std::string>(data);
    }

    /**
      * Copies the matrix from an existing \ref OptTaskSequence
      */
    OptTaskSequence& operator=(const OptTaskSequence& other);

    virtual ~OptTaskSequence(){
        tasks->clear();
        delete tasks;
    }

    /**
      * Adds a new task string to the sequence
      */
    void addTask(const std::string newtask) {tasks->push_back(newtask);}

    /**
      * Returns the tasks sequence
      */
    const std::vector<std::string>& getValue() const { return *tasks; }

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == TaskSequenceOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptTaskSequence
      */
    static OptTaskSequence* dynacast(GurlsOption* opt) {
        if (opt->isA(TaskSequenceOption) ){
            return static_cast<OptTaskSequence*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptTaskSequence
      */
    static const OptTaskSequence* dynacast(const GurlsOption* opt) {
        if (opt->isA(TaskSequenceOption) ){
            return static_cast<const OptTaskSequence*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Parse the task string at a given index and returns the task description and the task name
      */
    void getTaskAt(int index, std::string& taskdesc, std::string& taskname) {
        if (!isValid((*tasks)[index], taskdesc, taskname)){
            throw new gException(gurls::Exception_Invalid_TaskSequence);
        }
    }

    /**
      * Returns the number of tasks into the sequence
      */
    long int size(){
        return tasks->size();
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        int n = this->tasks->size();
        ar & n;
        for (int i = 0; i < n; i++){
            ar & tasks->at(i);
        }
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        int n = 0;
        ar & n;
        tasks->clear();
        for(int i = 0; i < n; i++){
            std::string s("");
            ar & s;
            this->tasks->push_back(s);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()


};

/**
  * Writes an OptMatrix to a stream
  */
template <typename T>
std::ostream& OptMatrix<T>::operator << (std::ostream& os){
    return os << std::endl << this->getValue();
}

}

#endif // _GURLS_OPT_H_
