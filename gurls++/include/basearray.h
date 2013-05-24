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

#ifndef _GURLS_BASEARRAY_H_
#define _GURLS_BASEARRAY_H_

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "exceptions.h"

namespace gurls {

/**
  * Maximum vector size for printing.
  * Verctor with a largest size will not be printed out
  */
const static int MAX_PRINTABLE_SIZE = 200;

/**
  *  \brief BaseArray is the base class for all classes implementing
  *  vectors and matrices as arrays of cells.
  *  \tparam T Cells type.
  */
template <typename T>
class BaseArray {

protected:

    T* data;                    ///< Pointer to the data buffer
    unsigned long size;         ///< Data buffer length
    bool isowner;               ///< Flag indicating whether vector has ownership of (and has to deallocate in destructor) the pointed buffer or not

    void alloc(unsigned long n); ///< Allocates the \c n elements data buffer

public:

    /**
      * Default constructor, initializes an empty vector
      */
    BaseArray() : data(0), size(0), isowner(false) { /* TO BE DISCUSSED*/}

    /**
      * Initializes a vector of n elements.
      */
    BaseArray(unsigned long n) {this->alloc(n);}

    /**
      * Initializes a vector copying elements from another vector
      */
    BaseArray(const BaseArray<T>& other);

    /**
      * Copies values form vector \c other to this one
      */
    BaseArray<T>& operator=(const BaseArray<T>& other);

    /**
      * Sets all elements of the vector to the value specified in \c val
      */
    BaseArray<T>& operator=(const T& val);

    /**
      * Destructor
      */
    ~BaseArray(){ if (this->isowner && data != NULL && size > 0) { delete[] this->data; } }

    /**
      * Copies \c n elements of a given vector \c v to this vector starting from \c start
      */
    void set(const T * v, unsigned long n, unsigned long start = 0);

    /**
      * Resizes the vector to length \c n
      */
    void resize(unsigned long n);

    /**
      * Copies \c n elements from this vector to a new buffer \c v
      */
    void asarray(T * v, unsigned long n) const;

    /**
      * Initializes the vector with pseudo-random values
      */
    void randomize();

    /**
      * Returns vector length
      */
    unsigned long getSize() const { return this->size; }

    /**
      * Returns a const pointer to the vector buffer
      */
    const T* getData() const { return this->data; }

    /**
      * Returns a non-const pointer to the vector buffer
      */
    T* getData() { return this->data; }

    /**
      * Returns a const pointer to the begin of vector buffer
      */
    const T* begin() const { return this->data(); }

    /**
      * Returns a non-const pointer to the begin of vector buffer
      */
    T* begin() { return this->data(); }

    /**
      * Returns a const pointer to the end of vector buffer
      */
    const T* end() const { return (this->data() + this->size); }

    /**
      * Returns a non-const pointer to the end of vector buffer
      */
    T* end() { return (this->data() + this->size); }


    /**
      * Returns a reference to the element with the largest value in the vector
      */
    const T& max() const;

    /**
      * Returns a reference to the element with the smallest value in the vector
      */
    const T& min() const;

    /**
      * Returns the sum of all elements in the vector
      */
    T sum() const;

    /**
      * Adds a value to all elements
      */
    BaseArray<T>& operator+=(T);

    /**
      * Subtracts a value to all elements
      */
    BaseArray<T>& operator-=(T);

    /**
      * Multiplies all elements by a scalar
      */
    BaseArray<T>& operator*=(T);

    /**
      * Divides all elements by a scalar
      */
    BaseArray<T>& operator/=(T);

    /**
      * In-place addition to a vector
      */
    BaseArray<T>& add(const BaseArray<T>&);

    /**
      * In-place subtraction to a vector
      */
    BaseArray<T>& subtract(const BaseArray<T>&);

    /**
      * In-place element by element multiplication by a vector
      */
    BaseArray<T>& multiply(const BaseArray<T>&);

    /**
      * In-place element by element division by a vector
      */
    BaseArray<T>& divide(const BaseArray<T>&);

    /**
      * In-place multiplicative inverse of each element
      */
    BaseArray<T>& setReciprocal();

    /**
      * Checks if all elements in a vector are equal to a given value
      */
    template <typename U>
    friend bool operator== (const BaseArray<U>&, const U&);

    /**
      * Checks if vector values are element by element close to a given array
      */
    bool closeTo(const BaseArray<T>&, T tolerance) const;

    /**
      * Returns a string description of the vector
      */
    virtual std::string what() const = 0;

    friend class boost::serialization::access;

    /**
      * Serializes the vector to a generic archive
      */
    template<class Archive>
    void save(Archive & , const unsigned int) const;

    /**
      * Deserializes the vector from a generic archive
      */
    template<class Archive>
    void load(Archive & , const unsigned int);

    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

}

#include "basearray.hpp"

#endif // _GURLS_BASEARRAY_H_

