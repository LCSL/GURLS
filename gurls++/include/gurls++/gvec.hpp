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

namespace gurls {


template <typename T>
gVec<T>::gVec(unsigned long n) {
    this->alloc(n);
}

template <typename T>
gVec<T>::gVec(T* buf, unsigned long n, bool owner) {
    if (owner) {
        this->alloc(n);
        this->set(buf, n);
    } else {
        this->isowner = owner;
        this->size = n;
        this->data = buf;
    };
}

template <typename T>
gVec<T>::gVec(const gVec<T>& other) {
    this->alloc(other.size);
    *this = other;
}

template <typename T>
gVec<T>& gVec<T>::operator=(const gVec<T>& other) {

    if (this != &other)
        this->set(other.data, other.size);

    return *this;
}

template <typename T>
gVec<T> gVec<T>::zeros(unsigned long n) {
    gVec<T> v(n);
    memset(v.data, 0, n*sizeof(T));
    return v;
}

template <typename T>
gVec<T> gVec<T>::rand(unsigned long n) {
    gVec<T> v(n);
    T* d = v.data;
    for (unsigned long i = 0; i < n; i++) {
        *d++ = (T)std::rand()/RAND_MAX;
    }
    return v;
}

template <typename T>
gVec<unsigned long> gVec<T>::nonzeros() const
{
    std::vector<unsigned long> indices;

    const T zero = static_cast<T>(0.0);

    for(unsigned long i=0; i<this->size; ++i)
    {
        if(!eq(this->data[i], zero))
            indices.push_back(i);
    }

    gVec<unsigned long> ret(indices.size());

    unsigned long *r_it = ret.getData();
    for(std::vector<unsigned long>::const_iterator it = indices.begin(), end = indices.end(); it != end; ++it, ++r_it)
        *r_it = *it;

    return ret;
}


template <typename T>
void gVec<T>::isequal(const T& value, std::vector<int>& indices) {
    const T* p = this->data;
    for (unsigned int i = 0; i < this->size; ++i, ++p) {
        if (*p == value)
            indices.push_back(i);
    }
}

/*
template <typename T>
void gVec<T>::resize(unsigned long n) {
  if (! this->m_owner) {
 assert(false);
  }
  else if (n == this->size) {
 return;
  } else {
 delete[] this->data;
 this->alloc(n);
  };
};
*/


template <typename T>
gVec<T> gVec<T>::subvec(unsigned int len, unsigned int start) const {
    gVec<T> v(len);
    v.set(&(this->data[start]), len);
    return v;
}


template <typename T>
gVec<T> gVec<T>::operator+(T val) const {
    gVec<T> w(*this);
    w += val;
    return w;
}

/**
  * Returns a vector containing the sum between a vector and a scalar
  */
template <typename T>
gVec<T> operator+(T val, const gVec<T>& v) {
    gVec<T> w(v);
    w += val;
    return w;
}

/**
  * Returns a vector containing the difference between a vector and a scalar
  */
template <typename T>
gVec<T> operator-(T val, const gVec<T>& v) {
    gVec<T> w(-v);
    w += val;
    return w;
}

// ------------ MULT SCALAR -------------------------------------------

/*
template <typename T>
gVec<T>& gVec<T>::operator*=(T val) {
  T* ptr = this->data;
  for (int i = 0; i < this->size; ++i, ++ptr) {
 *ptr *= val;
  }
  return *this;
};
*/

template <typename T>
gVec<T> gVec<T>::operator*(T val) const {
    gVec<T> w(*this);
    w *= val;
    return w;
}

/**
  * Returns a vector containing the multiplication of a vector by a scalar
  */
template <typename T>
gVec<T> operator*(T val, const gVec<T>& v) {
    gVec<T> w(v);
    w *= val;
    return w;
}

/*
template <typename T>
gVec<T>& gVec<T>::operator/=(T val) {
  std::cout << "bla" << std::endl;
  T* ptr = this->data;
  for (int i = 0; i < this->size; ++i, ++ptr) {
 *ptr /= val;
  }
  return *this;
};
*/

/**
  * Returns a vector containing the division of a vector by a scalar
  */
template <typename T>
gVec<T> operator/(T val, const gVec<T>& v) {
    gVec<T> w(v);
    w *= (static_cast<T>(1)/val);
    return w;
}

// ----------------- SUM OF VECTORS --------------------------

/*
template <typename T>
gVec<T>& gVec<T>::operator+=(const gVec<T>& v) {
  T *ptr = this->data, *ptr_v = v.data;
  for (int i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
 *ptr += *ptr_v;
  }
  return *this;
}
*/

template <typename T>
gVec<T>& gVec<T>::operator+=(const gVec<T>& v) {
    return static_cast<gVec<T>&>(BaseArray<T>::add(v));
}

template <typename T>
gVec<T> gVec<T>::operator+(const gVec<T>& v) const {
    gVec<T> w(v);
    w += *this;
    return w;
}

// ----------------- SUBTRACTION OF VECTORS --------------------------

template <typename T>
gVec<T>& gVec<T>::operator-=(const gVec<T>& v) {
    /*
  T *ptr = this->data, *ptr_v = v.data;
  for (int i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
 *ptr -= *ptr_v;
  }
  return *this;
  */
    return static_cast<gVec<T>&>(BaseArray<T>::subtract(v));
}

template <typename T>
gVec<T> gVec<T>::operator-(const gVec<T>& v) const {
    gVec<T> w(*this);
    w -= v;
    return w;
}

// ------------------- MULT TWO VECTORS ------------------------------------

template <typename T>
gVec<T>& gVec<T>::operator*=(const gVec<T>& v) {
    /*
  T *ptr = this->data, *ptr_v = v.data;
  for (int i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
 *ptr *= *ptr_v;
  }
  return *this;
  */
    return static_cast<gVec<T>&>(BaseArray<T>::multiply(v));
}

template <typename T>
gVec<T> gVec<T>::operator*(const gVec<T>& v) const {
    gVec<T> w(v);
    w *= *this;
    return w;
}

// ------------------- DIVIDE TWO VECTORS ------------------------------------

template <typename T>
gVec<T>& gVec<T>::operator/=(const gVec<T>& v) {
    /*
  T *ptr = this->data, *ptr_v = v.data;
  for (int i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
 *ptr /= *ptr_v;
  }
  return *this;
  */
    return static_cast<gVec<T>&>(BaseArray<T>::divide(v));
}

template <typename T>
gVec<T> gVec<T>::operator/(const gVec<T>& v) const {
    gVec<T> w(v);
    w /= *this;
    return w;
}


// TODO: to implement the comparison between two vectors!!

/**
  * Checks if all elements in a vector are equal to a given value
  */
template <typename U>
bool operator== (const gVec<U>& v, const U& val) {
    U *ptr = v.data;
    for (int i = 0; i < v.size; ++i, ++ptr) {
        if (*ptr != val){
            return false;
        }
    }
    return true;
}


// ----------------------- PRINTOUT --------------------------------

/**
  * Writes vector's information and data to a stream
  */
template <typename T>
std::ostream& operator<<(std::ostream& os, const gVec<T>& v) {
    if (v.getSize() >= (unsigned long)gurls::MAX_PRINTABLE_SIZE){
        os << v.what() << std::endl;
        return os;
    }
    os << "[" << v.data[0];
    for (unsigned long i = 1; i < v.size; ++i) {
        os << " " << v.data[i];
    }
    os << "]" << std::endl;
    return os;
}

template <typename T>
gVec<T>& gVec<T>::reciprocal() const {
    gVec<T> *w = new gVec(*this);
    w->setReciprocal();
    return *w;
}


template <typename T>
T gVec<T>::argmin() const {
    const T* val = std::min_element(this->data, this->data+this->size);
    return static_cast<T>(val-this->data);
}

template <typename T>
T gVec<T>::argmax() const {
    const T* val = std::max_element(this->data, this->data+this->size);
    return static_cast<T>(val-this->data);
}

template <typename T>
gVec<T>& gVec<T>::copyLocations(const gVec<T> locs){

    gVec<T>* v = new gVec<T>(locs.getSize());
    const T* ptr_locs = locs.getData();
    T* ptr_v = v->data;
    unsigned long val;
    for (unsigned long i = 0; i < locs.getSize(); i++) {
        val = static_cast<int>(*(ptr_locs+i));
        if ( (val < 0) || (val > this->size) ){
            throw gException(gurls::Exception_Index_Out_of_Bound);
        }
        *(ptr_v+i) = *(this->data+val);
    }
    return *v;
}

//template <typename T>
//template<class Archive>
//void gVec<T>::save(Archive & ar, const unsigned int /* file_version */) const{
//	ar & this->numrows & this->numcols & this->isowner;
//	T* ptr = this->data;
//	T* ptr_end = this->data+(this->numrows*this->numcols);
//	while (ptr!=ptr_end){
//		ar & *ptr++;
//	}
//}

//template <typename T>
//template<class Archive>
//void gMat2D<T>::load(Archive & ar, const unsigned int /* file_version */){
//	ar & this->numrows;
//	ar & this->numcols;
//	ar & this->isowner;
//	this->size = this->numrows*this->numcols;
//	this->data = new T[this->size];
//	T* ptr = this->data;
//	T* ptr_end = this->data+this->size;
//	while (ptr!=ptr_end){
//		ar & *ptr++;
//	}
//}



}
