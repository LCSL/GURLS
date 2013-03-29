/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, IIT@MIT Lab
 * All rights reserved.
 *
 * authors:  M. Santoro
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


#ifndef _GURLS_BIGARRAY_H_
#define _GURLS_BIGARRAY_H_

#include <hdf5.h>

#include <gmat2d.h>
#include <string>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/signal.hpp>
#include <boost/signals/connection.hpp>

namespace gurls
{


template<typename T>
hid_t getHdfType()
{
    throw gException(Exception_Unsupported_MatrixType);
}

template<>
GURLS_EXPORT hid_t getHdfType<float>();

template<>
GURLS_EXPORT hid_t getHdfType<double>();

template<>
GURLS_EXPORT hid_t getHdfType<unsigned long>();

template<>
GURLS_EXPORT hid_t getHdfType<unsigned int>();


template<typename T>
class BigArray: protected gMat2D<T>
{

public:

    BigArray();

    BigArray(std::string fileName);

    BigArray(std::string fileName, unsigned long r, unsigned long c);

    BigArray(std::string fileName, const gMat2D<T>& mat);

    BigArray(std::wstring fileName);

    BigArray(std::wstring fileName, unsigned long r, unsigned long c);

    BigArray(std::wstring fileName, const gMat2D<T>& mat);

    BigArray(const BigArray<T>& other);

//    BigArray(const BigArray<T>& other, std::string fileName);

//    BigArray<T>& operator= (BigArray<T> other);

    ~BigArray();

    void flush();

    void close();

    unsigned long rows() const;

    unsigned long cols() const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream&, const BigArray<U>&);

    virtual std::string what() const;


    // Getters
    gVec<T> operator[](unsigned long i) const;

    gVec<T> getRow(unsigned long i);

    gVec<T> operator() (unsigned long i) const;

    gVec<T> getColumn(unsigned long i) const;

    T operator() (unsigned long row, unsigned long col) const;

    T getValue(unsigned long row, unsigned long col) const;


    void getMatrix(unsigned long startingRow, unsigned long startingCol, unsigned long numRows, unsigned long numCols, gMat2D<T>&result) const;

    void getMatrix(unsigned long startingRow, unsigned long startingCol, gMat2D<T>&result) const;


    // Setter
    void setValue(unsigned long row, unsigned long col, T value);

    void setMatrix(unsigned long startingRow, unsigned long startingCol, const gMat2D<T>&value);

    void setMatrix(unsigned long startingRow, unsigned long startingCol, const T* M, const unsigned long M_rows, const unsigned long M_cols);

    void setColumn (unsigned long col, const gVec<T>& value);

    void setRow(unsigned long row, const gVec<T>& value);

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
      * Read the bigarray from a CSV file
      */
    void readCSV(const std::string& fileName);

    /**
      * Read the bigarray from a CSV file
      */
    void readCSV(const std::wstring& fileName);

    static void releaseMPIData()
    {
        BigArray<T>::releaseSignal()();
    }

    static bool isBigArray(const gMat2D<T>* matrix)
    {
        return (dynamic_cast<BigArray<T>* >(matrix) != NULL);
    }

    static bool isBigArray(const gMat2D<T>& matrix)
    {
        return isBigArray(&matrix);
    }

protected:

    static boost::signal0<void> &releaseSignal()
    {
        static boost::signal0<void> signal;
        return signal;
    }

    void loadNC(const std::string &fileName);

    void init(std::string& fileName, unsigned long r, unsigned long c);

    void getMatrix(unsigned long startingRow, unsigned long startingCol, unsigned long numRows, unsigned long numCols, T* result) const;


    hid_t file_id;
    hid_t dset_id;
    hid_t plist_id;

    std::string dataFileName;

    boost::signals::connection connection;
};

}

#include<bigarray.hpp>

#endif // _GURLS_BIGARRAY_H_
