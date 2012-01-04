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

#ifndef  _GURLS_GVEC_H_
#define  _GURLS_GVEC_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <cassert>
#include <queue>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "basearray.h"
#include "gmath.h"
#include "exceptions.h"


namespace gurls {

template <typename T>
class gMat2D;

// Implements a vector of generic length.
template <typename T>
class gVec : public BaseArray<T> {

public:
	gVec(unsigned long n = 0);
	gVec(T* buf, unsigned long n, bool ownership = true);
	gVec(const gVec& other);
	gVec<T>& operator=(const gVec& other);

	static gVec<T> zeros(unsigned long n = 0);

	static gVec<T> rand(unsigned long n = 0);

	gVec<unsigned long> nonzeros() const;

	// Make the vector as if it had been default-constructed.
	void clear();

	//using BaseArray<T>::size;

	gVec<T>& operator=(const T& val) {
		return static_cast<gVec<T>&>(this->BaseArray<T>::operator=(val));
	}

	//gVec<T>& set(const T * v, unsigned long n, unsigned long start = 0);
	using BaseArray<T>::set;

	//gVec<T>& set(const gVec<T>& v, unsigned long start=0) { return this->set(v.data(), v.size(), start); }
	void set(const gVec<T>& v, unsigned long start=0) { this->set(v.data(), v.size(), start); }

	//using BaseArray<T>::asarray;
	//void asarray(T * v, unsigned long n) const;

	const T& at(unsigned long i) const { return this->data[i]; }

	T& at(unsigned long i) { return this->data[i]; }

	const T& operator[](unsigned long i) const {
#ifdef DEBUG
		if (i < 0 || i >= this->size) {
			throw gException("Index exceeds array dimension.");
		}
#endif // DEBUG
		return at(i);
	};

	T& operator[](unsigned long i) {
#ifdef DEBUG
		if (i < 0 || i >= this->size) {
			throw gException("Index exceeds array dimension.");
		}
#endif // DEBUG
		return at(i);
	};

	gVec<T> subvec(unsigned int len, unsigned int start=0) const;

	gVec<T> operator-() const { return static_cast<T>(-1)*(*this); }

	using BaseArray<T>::operator+=;
	gVec<T> operator+(T) const;
	template <typename U>
	friend gVec<U> operator+(U val, const gVec<U>& v);

	using BaseArray<T>::operator-=;
	gVec<T> operator-(T val) const { return *this + (-val); }
	template <typename U>
	friend gVec<U> operator-(U val, const gVec<U>& v);

	//gVec<T>& operator*=(T );
	using BaseArray<T>::operator*=;
	gVec<T> operator*(T) const;
	template <typename U>
	friend gVec<U> operator*(U val, const gVec<U>& v);

	//gVec<T>& operator/=(T val);
	using BaseArray<T>::operator/=;
	gVec<T> operator/(T val) const { return *this * (static_cast<T>(1)/val); }
	template <typename U>
	friend gVec<U> operator/(U val, const gVec<U>& v);

	gVec<T>& operator+=(const gVec<T>&);
	gVec<T> operator+(const gVec<T>&) const;

	gVec<T>& operator-=(const gVec<T>& v);
	gVec<T> operator-(const gVec<T>& v) const;

	gVec<T>& operator*=(const gVec<T>&);
	gVec<T> operator*(const gVec<T>&) const;

	gVec<T>& operator/=(const gVec<T>&);
	gVec<T> operator/(const gVec<T>&) const;

	gVec<T>& reciprocal() const;


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

	void isequal(const T& value, std::vector<int>& indices);

	template <typename U>
	friend bool operator== (const gVec<U>&, const U&);

	template <typename U>
	friend std::ostream& operator<<(std::ostream&, const gVec<U>&);

	template <typename U>
	friend void dot(const gMat2D<U>&, const gVec<U>&, gVec<U>&);

	template <typename U>
	friend void dot(const gVec<T>& x, const gVec<T>& y);

	virtual std::string what() const {
		std::stringstream v("gVec:");
		v << std::string(" vector of length ");
		v << this->getSize();
		v <<  std::string("and type ");
		v << typeid(T).name();
		return v.str();
	}

	T argmin() const;
	T argmax() const;

	gVec<T>& copyLocations(const gVec<T> locs);



//	friend class boost::serialization::access;
//	template<class Archive>
//	void save(Archive & , const unsigned int) const;
//	template<class Archive>
//	void load(Archive & , const unsigned int);
//	BOOST_SERIALIZATION_SPLIT_MEMBER()

};


}

#include "gvec.hpp"

#endif //  _GURLS_GVEC_H_
