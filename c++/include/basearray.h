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

using namespace boost;

namespace gurls {

const static int MAX_PRINTABLE_SIZE = 200;

template <typename T>
class BaseArray {

protected:

	T* data;
	unsigned long size;
	bool isowner;

	void alloc(unsigned long n);

public:

	BaseArray() : data(0), size(0), isowner(false) { /* TO BE DISCUSSED*/};
	BaseArray(unsigned long n) {this->alloc(n);}
	BaseArray(const BaseArray& other);

	BaseArray<T>& operator=(const BaseArray& other);
	BaseArray<T>& operator=(const T& val);

	~BaseArray(){ if (this->isowner) { delete[] this->data; } };

	void set(const T * v, unsigned long n, unsigned long start = 0);
	void resize(unsigned long n);
	void asarray(T * v, unsigned long n) const;
	void randomize();

	unsigned long getSize() const { return this->size; }
	const T* getData() const { return this->data; }
	T* getData() { return this->data; }

	const T* begin() const { return this->data(); }
	T* begin() { return this->data(); }
	const T* end() const { return (this->data() + this->size); }
	T* end() { return (this->data() + this->size); }

	const T& max() const;
	const T& min() const;
	const double& sum() const;

	BaseArray<T>& operator+=(T);
	BaseArray<T>& operator-=(T);
	BaseArray<T>& operator*=(T);
	BaseArray<T>& operator/=(T);

	// The following methods implements in-place arithmetic operations between arrays
	BaseArray<T>& add(const BaseArray<T>&);
	BaseArray<T>& subtract(const BaseArray<T>&);
	BaseArray<T>& multiply(const BaseArray<T>&);
	BaseArray<T>& divide(const BaseArray<T>&);

	// In-place multiplicative inverse of each element
	BaseArray<T>& setReciprocal();

	template <typename U>
	friend bool operator== (const BaseArray<U>&, const U&);
	bool closeTo(const BaseArray<T>&, T tolerance) const;

	virtual std::string what() const = 0 ;

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & , const unsigned int) const;
	template<class Archive>
	void load(Archive & , const unsigned int);
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

}

#include "basearray.hpp"

#endif // _GURLS_BASEARRAY_H_

