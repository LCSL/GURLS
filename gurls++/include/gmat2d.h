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

#ifndef  _GURLS_GMAT2D_H_
#define  _GURLS_GMAT2D_H_


#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <queue>
#include <limits>
#include <typeinfo>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "basearray.h"
#include "exceptions.h"
#include "gvec.h"
#include "gmath.h"

namespace gurls {

static const int COLUMNWISE = 0; ///< Used to tell methods to operate on a matrix in column-wise order
static const int ROWWISE = 1;    ///< Used to tell methods to operate on a matrix in row-wise order


template <typename T>
class gVec;

/**
  * \ingroup LinearAlgebra
  * \brief gMat2D implements a matrix of generic size
  * \tparam T Cells type.
  */
template <typename T>
class gMat2D : public BaseArray<T> {

protected:

    unsigned long numcols;  ///< Number of columns
    unsigned long numrows;  ///< Number of rows

public:

    /**
      * Initializes a r-by-c matrix
      */
    gMat2D(unsigned long r = 0, unsigned long c = 0);

    /**
      * Initializes an r-by-c matrix from a data buffer. If ownership == true the class will take the ownership of the buffer
      */
    gMat2D(T* buf, unsigned long r, unsigned long c, bool owner);

    /**
      * Constructor form an existing \ref gMat2D
      */
    gMat2D(const gMat2D<T>& other);

    /**
      * Assignment operator
      */
    gMat2D<T>& operator=(const gMat2D<T>& other);

    /**
      * Returns a r-by-c matrix of all zeros
      */
    static gMat2D<T> zeros(unsigned long r = 0, unsigned long c = 0);

    /**
      * Returns a r-by-c matrix initialized with pseudo-random values
      */
    static gMat2D<T> rand(unsigned long r = 0, unsigned long c = 0);

    /**
      * Returns a squared matrix initialized in the diagonal with values from \c diagonal
      */
    static gMat2D<T> diag(gVec<T> diagonal);

    /**
      * Returns a n-by-n identity matrix
      */
    static gMat2D<T> eye(unsigned long n = 0);

    /**
      * Resizes the matrix to r-by-c size
      */
    void resize(unsigned long r, unsigned long c);

    /**
      * Reshapes the matrix to r-by-c size. Warning: r*c must be equals to \ref numrows * \ref numcols
      */
    void reshape(unsigned long r, unsigned long c);

    /**
      * Returns the number of columns
      */
    unsigned long cols() const {return numcols; }

    /**
      * Returns the number of rows
      */
    unsigned long rows() const {return numrows; }

    // The current implementation Uses the LU decomposition
    //T det();

    /**
      * Sets all elements of the matrix to the value specified in \c val
      */
    gMat2D<T>& operator=(const T& val) {
        return static_cast<gMat2D<T>&>(this->BaseArray<T>::operator=(val));
    }

    /**
      * Copies a r-by-c submatrix into a new matrix.
      * The input parameter \c other is supposed to be already initialized
      * with a number of rows equal to r and columns equal to c.
      */
    void submatrix(const gMat2D<T>& other, unsigned long r, unsigned long c);

    /**
      * Transpose the matrix.
      * The input parameter \c transposed is supposed to be already initialized
      * with a number of rows and columns equal to the number of columns and rows,
      * respectively, of this matrix.
      */
    void transpose(gMat2D <T>& transposed) const;

    /**
      * Returns a vector containing the linearized matrix
      */
    gVec<T> asvector() const {
        return gVec<T>(this->data, this->size, true);
    }

    /**
      * Returns a vector containing the elements at the i-th row
      */
    gVec<T> operator[](unsigned long i) const {
        return gVec<T>(this->data+i*this->numcols, this->numcols, true);
    }

    /**
      * Returns a vector containing the elements at the i-th row
      */
    gVec<T> operator[](unsigned long i) {
        return gVec<T>(this->data+i*this->numcols, this->numcols, false);
    }

    //-------------------------------------------------------------------------

    /**
      * Provides access to the elements of the matrix in a Matlab style: M(r,c)
      */
    T& operator() (unsigned long row, unsigned long col)
    {
        return this->data[this->cols()*row + col];
    }

    /**
      * Returns a vector containing the elements at the i-th column
      */
    gVec<T> operator() (unsigned long i)
    {
        gVec<T> v(this->numrows);
        for(unsigned long j=0;j<this->numrows;j++)
        {
            v[j]=this->data[this->cols()*j+i];
        }
        return v;
    }

    /**
      * Returns a vector containing the elements at the i-th column
      */
    gVec<T> operator() (unsigned long i) const
    {
        gVec<T> v(this->numrows);
        for(unsigned long j=0;j<this->numrows;j++)
        {
            v[j]=this->data[this->cols()*j+i];
        }
        return v;
    }

    /**
      * Sets the elements at the i-th row
      */
    void setRow(gVec<T>& vec, unsigned long i){
        T* ptr = this->data+this->numcols*i;
        T* ptr_end = this->data+this->numcols*(i+1);
        T* v_ptr = vec.getData();
        while (ptr != ptr_end){
            *ptr++ = *v_ptr++;
        }
        return;
    }

    /**
      * Sets the elements at the i-th columns
      */
    void setColumn(gVec<T>& vec, unsigned long i){
        for(int j=0;j<this->numrows;j++)
        {
            this->data[this->cols()*j+i]=vec[j];
        }
    }

    /**
      * Sets the elements on the matrix diagonal
      */
    void setDiag(gVec<T>& vec){
        unsigned long int k = std::min(this->numcols, this->numrows);
        if (vec.getSize() < k) {
            throw gException("The lenght of the vector must be at least equal to the minimum dimension of the matrix");
        }
        T* ptr = this->data;
        T* ptrVec = vec.getData();

        for(unsigned long int j = 0; j < k; j++, ptrVec++, ptr+=(this->numcols+1)){
            *ptr = *ptrVec;
        }
    }

    /**
      * Inverts elements sign
      */
    gMat2D<T> operator-() const { return static_cast<T>(-1)*(*this); }

    using BaseArray<T>::operator+=;

    /**
      * Returns a matrix containing the sum between the matrix and a scalar
      */
    gMat2D<T> operator+(T) const;

    /**
      * Returns a matrix containing the sum between a matrix and a scalar
      */
    template <typename U>
    friend gMat2D<U> operator+(U val, const gMat2D<U>& v);

    using BaseArray<T>::operator-=;

    /**
      * Returns a matrixr containing the difference between the matrix and a scalar
      */
    gMat2D<T> operator-(T val) const { return *this + (-val); }

    /**
      * Returns a matrix containing the difference between a matrix and a scalar
      */
    template <typename U>
    friend gMat2D<U> operator-(U val, const gMat2D<U>& v);


    using BaseArray<T>::operator*=;
    /**
      * Returns a matrix containing the multiplication of the matrix by a scalar
      */
    gMat2D<T> operator*(T) const;

    /**
      * Returns a matrix containing the multiplication of a matrix by a scalar
      */
    template <typename U>
    friend gMat2D<U> operator*(U val, const gMat2D<U>& v);

    using BaseArray<T>::operator/=;

    /**
      * Returns a matrix containing the division of the matrix by a scalar
      */
    gMat2D<T> operator/(T val) const { return *this * (static_cast<T>(1)/val); }

    /**
      * Returns a matrix containing the division of a matrix by a scalar
      */
    template <typename U>
    friend gMat2D<U> operator/(U val, const gMat2D<U>& v);

    /**
      * In-place element-by-element addition to a matrix
      */
    gMat2D<T>& operator+=(const gMat2D<T>&);

    /**
      * Returns the element-by-element sum of two vectors
      */
    gMat2D<T> operator+(const gMat2D<T>&) const;

    /**
      * In-place element-by-element subtraction to a matrix
      */
    gMat2D<T>& operator-=(const gMat2D<T>& v);

    /**
      * Returns the element-by-element difference between two matrices
      */
    gMat2D<T> operator-(const gMat2D<T>& v) const;

    /**
      * In-place element-by-element multiplication by a matrix
      */
    gMat2D<T>& operator*=(const gMat2D<T>&);

    /**
      * Returns the element-by-element multiplication between two matrices
      */
    gMat2D<T> operator*(const gMat2D<T>&) const;

    /**
      * In-place element-by-element division by a matrix
      */
    gMat2D<T>& operator/=(const gMat2D<T>&);

    /**
      * Returns the element-by-element division between two matrices
      */
    gMat2D<T> operator/(const gMat2D<T>&) const;

    /**
      * Returns true if all elements of the matrix are equals to a value \c val
      */
    bool allEqualsTo(const T& val)  const;

    /**
      * Writes matrix information and data to a stream
      */
    template <typename U>
    friend std::ostream& operator<<(std::ostream&, const gMat2D<U>&);

    /**
      * General Matrix-Matrix multiplication of two matrices
      */
    template <typename U>
    friend void dot(const gMat2D<U>&, const gMat2D<U>&, gMat2D<U>&);

    /**
      * General Matrix-Vector multiplication of a matrix with a vector
      */
    template <typename U>
    friend void dot(const gMat2D<U>&, const gVec<U>&, gVec<U>&);

    // ================= TO BE DISCUSSED =======================
//    template <typename U>
//    gMat2D<U> repmat(const gVec<U>&, unsigned long, bool);
    // =========================================================

    /**
      * Clears the lower triangle of a matrix
      */
    void uppertriangular(gMat2D<T>& up) const;

    /**
      * Clears the upper triangle of a matrix
      */
    void lowertriangular(gMat2D<T>& lo) const;


    static const std::string Less;          ///< String for < operator
    static const std::string Greater;       ///< String for > operator
    static const std::string LessEq;        ///< String for <= operator
    static const std::string GreaterEq;     ///< String for >= operator
    static const std::string Equal;         ///< String for == operator

    /**
      * Compares each element of the matrix with the threshold
      * and returns a matrix with elements set to true where element == threshold,
      * or false where element != threshold
      */
    gMat2D<bool> & operator ==(T threshold) const;

    /**
      * Compares each element of the matrix with the threshold
      * and returns a matrix with elements set to true where element < threshold,
      * or false where element >= threshold
      */
    gMat2D<bool> & operator <(T threshold) const;

    /**
      * Compares each element of the matrix with the threshold
      * and returns a matrix with elements set to true where element <= threshold,
      * or false where element > threshold
      */
    gMat2D<bool> & operator <=(T threshold) const;

    /**
      * Compares each element of the matrix with the threshold
      * and returns a matrix with elements set to true where element > threshold,
      * or false where element <= threshold
      */
    gMat2D<bool> &  operator >(T threshold) const;

    /**
      * Compares each element of the matrix with the threshold
      * and returns a matrix with elements set to true where element >= threshold,
      * or false where element < threshold
      */
    gMat2D<bool> & operator >=(T threshold) const;

    /**
      * Compares each element of the matrix with the threshold
      * and returns a matrix with elements set to true or false depending on the
      * comparison operator passed as parameter
      */
    gMat2D<bool> & compare(T& threshold, std::string logical_operator = GreaterEq) const;

    /**
      * Returns a vector containing all matrix elements where comparison matrix
      * == true
      */
    gVec<T>& where(const gMat2D<bool> & comparison) const;

    /**
      * Returns matrix element-by-element multiplicative inverse
      */
    gMat2D<T>& reciprocal() const;

    /**
      * Returns a string description of the matrix
      */
    virtual std::string what() const {
        std::stringstream v("gMat2D:");
        v << std::string(" (");
        v << this->rows() ;
        v << std::string(" x ");
        v << this->cols();
        v <<std::string(") matrix of type ");
        v << typeid(T).name();
        return v.str();
    }

    using BaseArray<T>::max;
    using BaseArray<T>::min;
    using BaseArray<T>::sum;

    /**
      * Returns a vector containing the smallest elements along the columns
      * if order == \ref COLUMNWISE  or along the rows if order == \ref ROWWISE
      */
    gVec<T>& min(int order);

    /**
      * Returns a vector containing the largest elements along the columns
      * if order == \ref COLUMNWISE  or along the rows if order == \ref ROWWISE
      */
    gVec<T>& max(int order);

    /**
      * Returns a vector containing the smallest elements along the columns
      * if order == \ref COLUMNWISE  or along the rows if order == \ref ROWWISE
      */
    gVec<T>& argmin(int order);

    /**
      * Returns a vector containing the largest elements along the columns
      * if order == \ref COLUMNWISE  or along the rows if order == \ref ROWWISE
      */
    gVec<T>& argmax(int order);

    /**
      * Returns a vector containing the sums of the elements along the columns
      * if order == \ref COLUMNWISE  or along the rows if order == \ref ROWWISE
      */
    gVec<T>& sum(int order) const;

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

    /**
      * Serializes the matrix to file
      */
    void save(const std::string& fileName) const;

    /**
      * Deserializes the matrix from file
      */
    void load(const std::string& fileName);

};

}

#include "gmat2d.hpp"


#endif // _GURLS_GMAT2D_H_
