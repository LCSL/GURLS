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

static const int COLUMNWISE = 0;
static const int ROWWISE = 1;


template <typename T>
class gVec;

template <typename T>
class gMat2D : public BaseArray<T> {

protected:

	unsigned long numcols;
	unsigned long numrows;


public:

	gMat2D(unsigned long r = 0, unsigned long c = 0);
	gMat2D(T* buf, unsigned long r, unsigned long c, bool owner);
	gMat2D(const gMat2D& other);
	gMat2D<T>& operator=(const gMat2D& other);

	static gMat2D<T> zeros(unsigned long r = 0, unsigned long c = 0);

	static gMat2D<T> rand(unsigned long r = 0, unsigned long c = 0);

	static gMat2D<T> diag(gVec<T> diagonal);

	static gMat2D<T> eye(unsigned long n = 0);

	void resize(unsigned long r, unsigned long c);
	void reshape(unsigned long r, unsigned long c);

	unsigned long cols() const {return numcols; }

	unsigned long rows() const {return numrows; }

	// The current implementation Uses the LU decomposition
	T det();

	gMat2D<T>& operator=(const T& val) {
		return static_cast<gMat2D<T>&>(this->BaseArray<T>::operator=(val));
	}

	void submatrix(const gMat2D<T>& other, unsigned long r, unsigned long c);

	// The input argument ``transposed'' is supposed to be already initialized
	// with a number of rows and columns equal to the number of columns and rows,
	// respectively, of this matrix.
	void transpose(gMat2D <T>& transposed) const;

	gVec<T> asvector() const {
		return gVec<T>(this->data, this->size, true);
	};

	gVec<T> operator[](unsigned long i) const {
		return gVec<T>(this->data+i*this->numcols, this->numcols, true);
	};

	gVec<T> operator[](unsigned long i) {
		return gVec<T>(this->data+i*this->numcols, this->numcols, false);
	};

	//-------------------------------------------------------------------------
	// provides access to the elements of the matrix in a Matlab style: M(r,c)
	T& operator() (unsigned long row, unsigned long col)
	{
		return this->data[this->cols()*row + col];
	}

	// provides access to the i-th column of the matrix
	gVec<T> operator() (unsigned long col)
	{
		gVec<T> v(this->numrows);
		for(int j=0;j<this->numrows;j++)
		{
			v[j]=this->data[this->cols()*j+col];
		}
		return v;
	}

	void setRow(gVec<T>& vec, unsigned long row){
		T* ptr = this->data+this->numcols*row;
		T* ptr_end = this->data+this->numcols*(row+1);
		T* v_ptr = vec.getData();
		while (ptr != ptr_end){
			*ptr++ = *v_ptr++;
		}
		return;
	}

	void setColumn(gVec<T>& vec, unsigned long col){
		for(int j=0;j<this->numrows;j++)
		{
			this->data[this->cols()*j+col]=vec[j];
		}
	}

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

	gMat2D<T> operator-() const { return static_cast<T>(-1)*(*this); }

	//gMat2D<T>& operator+=(T );
	using BaseArray<T>::operator+=;
	gMat2D<T> operator+(T) const;
	template <typename U>
	friend gMat2D<U> operator+(U val, const gMat2D<U>& v);

	//gMat2D<T>& operator-=(T val) { return *this += (-val); }
	using BaseArray<T>::operator-=;
	gMat2D<T> operator-(T val) const { return *this + (-val); }
	template <typename U>
	friend gMat2D<U> operator-(U val, const gMat2D<U>& v);

	//gMat2D<T>& operator*=(T );
	using BaseArray<T>::operator*=;
	gMat2D<T> operator*(T) const;
	template <typename U>
	friend gMat2D<U> operator*(U val, const gMat2D<U>& v);

	//gMat2D<T>& operator/=(T val);
	using BaseArray<T>::operator/=;
	gMat2D<T> operator/(T val) const { return *this * (static_cast<T>(1)/val); }
	template <typename U>
	friend gMat2D<U> operator/(U val, const gMat2D<U>& v);

	gMat2D<T>& operator+=(const gMat2D<T>&);
	gMat2D<T> operator+(const gMat2D<T>&) const;

	gMat2D<T>& operator-=(const gMat2D<T>& v);
	gMat2D<T> operator-(const gMat2D<T>& v) const;

	gMat2D<T>& operator*=(const gMat2D<T>&);
	gMat2D<T> operator*(const gMat2D<T>&) const;

	gMat2D<T>& operator/=(const gMat2D<T>&);
	gMat2D<T> operator/(const gMat2D<T>&) const;

	bool allEqualsTo(const T& val)  const;

	template <typename U>
	friend std::ostream& operator<<(std::ostream&, const gMat2D<U>&);

	template <typename U>
	friend void dot(const gMat2D<U>&, const gMat2D<U>&, gMat2D<U>&);

	template <typename U>
	friend void dot(const gMat2D<U>&, const gVec<U>&, gVec<U>&);

	// ================= TO BE DISCUSSED =======================
	template <typename U>
	gMat2D<U> repmat(const gVec<U>&, unsigned long, bool);
	// =========================================================

	void uppertriangular(gMat2D<T>& up) const ;
	void lowertriangular(gMat2D<T>& lo) const ;

	static const std::string Less;
	static const std::string Greater;
	static const std::string LessEq;
	static const std::string GreaterEq;
	static const std::string Equal;

	gMat2D<bool> & operator ==(T threshold) const ;
	gMat2D<bool> & operator <(T threshold) const ;
	gMat2D<bool> & operator <=(T threshold) const ;
	gMat2D<bool> &  operator >(T threshold) const ;
	gMat2D<bool> & operator >=(T threshold) const ;
	gMat2D<bool> & compare(T& threshold, std::string logical_operator = GreaterEq) const ;
	gVec<T>& where(const gMat2D<bool> & comparison) const ;

	gMat2D<T>& reciprocal() const;

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

	gVec<T>& min(int order) ;
	gVec<T>& max(int order) ;
	gVec<T>& argmin(int order) ;
	gVec<T>& argmax(int order) ;
	gVec<T>& sum(int order) ;

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & , const unsigned int) const;
	template<class Archive>
	void load(Archive & , const unsigned int);
	BOOST_SERIALIZATION_SPLIT_MEMBER()

};

}

#include "gmat2d.hpp"


#endif // _GURLS_GMAT2D_H_
