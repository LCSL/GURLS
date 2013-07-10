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

#ifndef _GURLS_OPTFUNCTION_H_
#define _GURLS_OPTFUNCTION_H_

#include "gurls++/exports.h"
#include "gurls++/options.h"
#include "gurls++/exceptions.h"

#include <string>

namespace gurls
{

class GURLS_EXPORT Functor
{
public:
    virtual double operator()(const double *v, int n) = 0;
    virtual float operator()(const float *v, int n) = 0;
};

/**
 * Computes the mean value of a vector \c v of lenght \c n
 */
class GURLS_EXPORT Mean : public Functor
{
public:

    double operator()(const double *v, int n)
    {
        return mean<double>(v, n);
    }

    float operator()(const float *v, int n)
    {
        return mean<float>(v, n);
    }

protected:
    template <typename T>
    T mean(const T *v, int n)
    {
        T m = (T)0.0;

        for (int i = 0; i < n; ++i, ++v)
            m += *v;

        return m/n;
    }
};

/**
 * Computes the smallest element in a vector \c v of lenght \c n
 */
class GURLS_EXPORT Min : public Functor
{
public:
    double operator()(const double *v, int n)
    {
        return *std::min_element(v,v+n);
    }

    float operator()(const float *v, int n)
    {
        return *std::min_element(v,v+n);
    }
};


/**
 * Computes the largest element in a vector \c v of lenght \c n
 */
class GURLS_EXPORT Max : public Functor
{
public:
    double operator()(const double *v, int n)
    {
        return *std::max_element(v,v+n);
    }

    float operator()(const float *v, int n)
    {
        return *std::max_element(v,v+n);
    }
};

/**
 * Computes the median value of a vector \c v of lenght \c n
 */
class GURLS_EXPORT Median : public Functor
{
public:
    double operator()(const double *v, int n)
    {
        return median(v,n);
    }

    float operator()(const float *v, int n)
    {
        return median(v,n);
    }

protected:
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
};




#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4251)
#endif


/**
  * \ingroup Settings
  * \brief OptFunction is an option representing a pointer to a generic function
  * T (*function)(T* , int) operating over an array of floating point numbers.
  */
class GURLS_EXPORT OptFunction: public GurlsOption
{
protected:
    std::string name; ///< Function name
    Functor *f;       ///< Pointer to the functor containing the function to be called

public:

    /**
      * Empty constructor
      */
    OptFunction();

    /**
      * Constructor from a fuction name
      */
    OptFunction(std::string func_name);

    /**
      * Destructor
      */
    ~OptFunction();

    /**
      * Copies the option values from an existing \ref OptFunction
      */
    OptFunction& operator=(const OptFunction& other);

    /**
      * Returns the function name
      */
    std::string getName() const;

    /**
      * Executes the function over a buffer of length n, returning the result
      */
    template<typename T>
    T getValue(T* array, int n) const
    {
        return (*f)(array, n);
    }

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const;

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptFunction
      */
    static OptFunction* dynacast(GurlsOption* opt);

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptFunction
      */
    static const OptFunction* dynacast(const GurlsOption* opt);

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);


protected:

    /**
      * Sets the option value from a string representing a supported function name
      */
    void setValue(std::string func_name);

};


#ifdef _WIN32
#pragma warning(pop)
#endif


}

#endif //_GURLS_OPTFUNCTION_H_
