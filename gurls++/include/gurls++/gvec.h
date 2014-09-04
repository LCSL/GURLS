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

#ifndef  _GURLS_GVEC_H_
#define  _GURLS_GVEC_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <cassert>
#include <sstream>
#include <typeinfo>


#include "gurls++/basearray.h"
#include "gurls++/gmath.h"
#include "gurls++/exceptions.h"


namespace gurls {

template <typename T>
class gMat2D;

/**
  * \ingroup LinearAlgebra
  * \brief gVec implements a vector of generic length
  * \tparam T Cells type.
  */
template <typename T>
class gVec : public BaseArray<T> {

public:
    /**
      * Initializes a vector of n elements.
      */
    gVec(unsigned long n = 0);

    /**
      * Initializes a vector of n elements from a data buffer. If ownership == true the class will take the ownership of the buffer
      */
    gVec(T* buf, unsigned long n, bool ownership = true);

    /**
      * Constructor form an existing \ref gVec
      */
    gVec(const gVec<T>& other);

    /**
      * Assignment operator
      */
    gVec<T>& operator=(const gVec<T>& other);

    /**
      * Returns a vector of all zeros
      */
    static gVec<T> zeros(unsigned long n = 0);

    /**
      * Returns a n-length vector initialized with pseudo-random values
      */
    static gVec<T> rand(unsigned long n = 0);

    /**
      * Returns a subvector containing all nonzero values
      */
    gVec<unsigned long> nonzeros() const;

    /**
      * Clears the vector, making it as if it had been default-constructed.
      */
    void clear();

    /**
      * Sets all elements of the vector to the value specified in \c val
      */
    gVec<T>& operator=(const T& val) { return static_cast<gVec<T>&>(this->BaseArray<T>::operator=(val)); }

    using BaseArray<T>::set;

    /**
      * Copies all elements of a given vector \c v to the vector starting from \c start
      */
    void set(const gVec<T>& v, unsigned long start=0) { this->set(v.data, v.size, start); }

    /**
      * Returns a const reference to the i-th element
      */
    const T& at(unsigned long i) const { return this->data[i]; }

    /**
      * Returns a reference to the i-th element
      */
    T& at(unsigned long i) { return this->data[i]; }

    /**
      * Returns a const reference to the i-th element
      */
    const T& operator[](unsigned long i) const {
#ifdef DEBUG
        if (i < 0 || i >= this->size)
            throw gException("Index exceeds array dimension.");
#endif // DEBUG
        return at(i);
    }

    /**
      * Returns a reference to the i-th element
      */
    T& operator[](unsigned long i) {
#ifdef DEBUG
        if (i < 0 || i >= this->size)
            throw gException("Index exceeds array dimension.");
#endif // DEBUG
        return at(i);
    }

    /**
      * Returns a subvector containing all values beginning from \c start
      */
    gVec<T> subvec(unsigned int len, unsigned int start=0) const;

    /**
      * Inverts elements sign
      */
    gVec<T> operator-() const { return static_cast<T>(-1)*(*this); }

    using BaseArray<T>::operator+=;

    /**
      * Returns a vector containing the sum between the vector and a scalar
      */
    gVec<T> operator+(T) const;

    /**
      * Returns a vector containing the sum between a vector and a scalar
      */
    template <typename U>
    friend gVec<U> operator+(U val, const gVec<U>& v);


    using BaseArray<T>::operator-=;

    /**
      * Returns a vector containing the difference between the vector and a scalar
      */
    gVec<T> operator-(T val) const { return *this + (-val); }

    /**
      * Returns a vector containing the difference between a vector and a scalar
      */
    template <typename U>
    friend gVec<U> operator-(U val, const gVec<U>& v);

    using BaseArray<T>::operator*=;

    /**
      * Returns a vector containing the multiplication of the vector by a scalar
      */
    gVec<T> operator*(T) const;

    /**
      * Returns a vector containing the multiplication of a vector by a scalar
      */
    template <typename U>
    friend gVec<U> operator*(U val, const gVec<U>& v);

    using BaseArray<T>::operator/=;

    /**
      * Returns a vector containing the division of the vector by a scalar
      */
    gVec<T> operator/(T val) const { return *this * (static_cast<T>(1)/val); }

    /**
      * Returns a vector containing the division of a vector by a scalar
      */
    template <typename U>
    friend gVec<U> operator/(U val, const gVec<U>& v);

    /**
      * In-place addition to a vector
      */
    gVec<T>& operator+=(const gVec<T>&);

    /**
      * Returns the sum of two vectors
      */
    gVec<T> operator+(const gVec<T>&) const;

    /**
      * In-place subtraction to a vector
      */
    gVec<T>& operator-=(const gVec<T>& v);

    /**
      * Returns the difference between two vectors
      */
    gVec<T> operator-(const gVec<T>& v) const;

    /**
      * In-place element by element multiplication by a vector
      */
    gVec<T>& operator*=(const gVec<T>&);

    /**
      * Returns the element by element multiplication between two vectors
      */
    gVec<T> operator*(const gVec<T>&) const;

    /**
      * In-place element by element division by a vector
      */
    gVec<T>& operator/=(const gVec<T>&);

    /**
      * Returns the element by element division between two vectors
      */
    gVec<T> operator/(const gVec<T>&) const;

    /**
      * Returns vector's multiplicative inverse
      */
    gVec<T>& reciprocal() const;

    /**
      * Returns a matrix representing the vector
      */
    gMat2D<T>& asMatrix(bool ascolumn = true){
        unsigned long int r = this->size;
        unsigned long int c = 1;
        if (!ascolumn){
            r = 1;
            c = this->size;
        }
        gMat2D<T>* mat = new gMat2D<T>(this->data, r, c, true);
        return *mat;
    }

    /**
      * Writes into the vector \c indices the indices of the elements whose value is equal to the given one
      */
    void isequal(const T& value, std::vector<int>& indices);

    /**
      * Checks if all elements in a vector are equal to a given value
      */
    template <typename U>
    friend bool operator== (const gVec<U>&, const U&);

    /**
      * Writes vector's information and data to a stream
      */
    template <typename U>
    friend std::ostream& operator<<(std::ostream&, const gVec<U>&);

//    template <typename U>
//    friend void dot(const gMat2D<U>&, const gVec<U>&, gVec<U>&);

//    template <typename U>
//    friend void dot(const gVec<T>& x, const gVec<T>& y);

    /**
      * Returns a string description of the vector
      */
    virtual std::string what() const {
        std::stringstream v("gVec:");
        v << std::string(" vector of length ");
        v << this->getSize();
        v <<  std::string("and type ");
        v << typeid(T).name();
        return v.str();
    }

    /**
     * Returns the smallest element in the vector
     */
    T argmin() const;

    /**
     * Returns the largest element in the vector
     */
    T argmax() const;

    /**
     * Returns a subvector containing elements whose index in the vector is contained in \c locs
     */
    gVec<T>& copyLocations(const gVec<T> locs);

};


}

#include "gvec.hpp"

#endif //  _GURLS_GVEC_H_
