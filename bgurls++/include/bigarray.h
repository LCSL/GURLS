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

#include <netcdf>

#include <gmat2d.h>
#include <string>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

namespace gurls
{

template<typename T>
netCDF::NcType getNcType()
{
    throw gException(Exception_Unsupported_MatrixType);
}

template<>
netCDF::NcType getNcType<float>();

template<>
netCDF::NcType getNcType<double>();

template<>
netCDF::NcType getNcType<unsigned long>();

template<typename T>
class BigArray: protected gMat2D<T>
{
protected:
    typedef typename netCDF::NcFile NcFile;
    typedef typename netCDF::NcVar NcVar;
    typedef typename netCDF::NcDim NcDim;
    typedef typename netCDF::exceptions::NcException NcException;

public:

    BigArray(): gMat2D<T>(), dataFile(NULL), dataFileName("")
    {
    }

    BigArray(std::string fileName)
    {
        loadNC(fileName);
    }

    BigArray(std::string fileName, unsigned long r, unsigned long c)
    {
        init(fileName, r, c);
    }

    BigArray(std::string fileName, const gMat2D<T>& mat)
    {
        init(fileName, mat.rows(), mat.cols());
        setMatrix(0, 0, mat);
    }

    BigArray(const BigArray<T>& other)
    {
        loadNC(other.dataFileName);
    }

//    BigArray(const BigArray<T>& other, std::string fileName);

    ~BigArray()
    {
        delete dataFile;
    }


    unsigned long rows() const
    {
        return this->numrows;
    }

    unsigned long cols() const
    {
        return this->numcols;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream&, const BigArray<U>&);

    virtual std::string what() const
    {
        std::stringstream v;
        v << "BigArray: (" << this->numrows;
        v << " x " << this->numcols;
        v << ") matrix of type " << typeid(T).name();
        return v.str();
    }

    // Getter

    gVec<T> operator[](unsigned long i) const
    {
        if(i >= this->numrows)
            throw gException(Exception_Index_Out_of_Bound);

        gVec<T> ret(this->numcols);

        std::vector<size_t> start(2), count(2);

        start[0] = 0;
        start[1] = i;

        count[0] = this->numcols;
        count[1] = 1;

        matrix.getVar(start, count, ret.getData());

        return ret;
    }

    gVec<T> getRow(unsigned long i)
    {
        return (*this)[i];
    }

    gVec<T> operator() (unsigned long i) const
    {
        if(i >= this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        gVec<T> ret(this->numrows);

        std::vector<size_t> start(2), count(2);

        start[0] = i;
        start[1] = 0;

        count[0] = 1;
        count[1] = this->numrows;

        matrix.getVar(start, count, ret.getData());

        return ret;
    }

    gVec<T> getColumnn(unsigned long i)
    {
        return (*this)(i);
    }

    T operator() (unsigned long row, unsigned long col) const
    {
        if(row >= this->numrows || col >= this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        T ret;

        std::vector<size_t> start(2), count(2, 1);

        start[0] = col;
        start[1] = row;

        matrix.getVar(start, count, &ret);

        return ret;
    }

    T getValue(unsigned long row, unsigned long col) const
    {
        return (*this)(row, col);
    }


    void getMatrix(unsigned long startingRow, unsigned long startingCol, unsigned long numRows, unsigned long numCols, gMat2D<T>&result) const
    {
        if(startingRow >= this->numrows || startingCol >= this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        if(startingRow+numRows > this->numrows || startingCol+numCols > this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        result.resize(numRows, numCols);

        std::vector<size_t> start(2), count(2);

        start[0] = startingCol;
        start[1] = startingRow;

        count[0] = numCols;
        count[1] = numRows;

        matrix.getVar(start, count, result.getData());

    }

    void getMatrix(unsigned long startingRow, unsigned long startingCol, gMat2D<T>&result) const
    {
        if(startingRow >= this->numrows || startingCol >= this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        if(startingRow+result.rows() > this->numrows || startingCol+result.cols() > this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        std::vector<size_t> start(2), count(2);

        start[0] = startingCol;
        start[1] = startingRow;

        count[0] = result.cols();
        count[1] = result.rows();

        matrix.getVar(start, count, result.getData());

    }


    // Setter

    void setValue(unsigned long row, unsigned long col, T value)
    {
        if(row >= this->numrows || col >= this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        std::vector<size_t> start(2), count(2, 1);

        start[0] = col;
        start[1] = row;

        matrix.putVar(start, count, &value);
    }

    void setMatrix(unsigned long startingRow, unsigned long startingCol, const gMat2D<T>&value)
    {
        if(startingRow >= this->numrows || startingCol >= this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        if(startingRow+value.rows() > this->numrows || startingCol+value.cols() > this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        std::vector<size_t> start(2), count(2);

        start[0] = startingCol;
        start[1] = startingRow;

        count[0] = value.cols();
        count[1] = value.rows();

        matrix.putVar(start, count, value.getData());
    }

    void setColumn (unsigned long col, const gVec<T>& value)
    {
        if(col >= this->numcols)
            throw gException(Exception_Index_Out_of_Bound);

        if(value.getSize() != this->numrows)
            throw gException(Exception_Inconsistent_Size);

        std::vector<size_t> start(2), count(2);

        start[0] = col;
        start[1] = 0;

        count[0] = 1;
        count[1] = this->numrows;

        matrix.putVar(start, count, value.getData());
    }

    void setRow(unsigned long row, const gVec<T>& value)
    {
        if(row >= this->numrows)
            throw gException(Exception_Index_Out_of_Bound);

        if(value.getSize() != this->numcols)
            throw gException(Exception_Inconsistent_Size);

        std::vector<size_t> start(2), count(2);

        start[0] = 0;
        start[1] = row;

        count[0] = this->numcols;
        count[1] = 1;

        matrix.putVar(start, count, value.getData());
    }

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

protected:

    void loadNC(const std::string &fileName)
    {
        dataFileName = fileName;
        try
        {
            dataFile = new NcFile(fileName, NcFile::write/*, NcFile::nc4*/);

            this->numrows = dataFile->getDim("cols").getSize();
            this->numcols = dataFile->getDim("rows").getSize();

            matrix = dataFile->getVar("mat");
        }
        catch(NcException& e)
        {
            throw gException("Error opening file " + fileName + ": " + e.what());
        }
    }

    void init(std::string& fileName, unsigned long r, unsigned long c)
    {
        dataFileName = fileName;

        try
        {
            dataFile = new NcFile(fileName, NcFile::replace/*, NcFile::nc4*/);

            std::vector<NcDim> dims(2);

            // gurls++ uses column-major order while netCDF uses row-major
            dims[1] = dataFile->addDim("cols", r);
            this->numrows = r;

            dims[0] = dataFile->addDim("rows", c);
            this->numcols = c;

            matrix = dataFile->addVar("mat", getNcType<T>(), dims);

        }
        catch(NcException& e)
        {
            throw gException("Error opening file " + fileName + ": " + e.what());
        }
    }

    NcFile* dataFile;
    NcVar matrix;

    std::string dataFileName;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const BigArray<T>& m)
{
   if (m.rows() >= (unsigned long)gurls::MAX_PRINTABLE_SIZE || m.cols() >= (unsigned long) gurls::MAX_PRINTABLE_SIZE)
   {
       os << m.what() << std::endl;
       return os;
   }

    os << "[";
    for (unsigned long i = 0; i < m.numrows; ++i)
    {
        for (unsigned long j = 0; j < m.numcols; ++j)
            os << " " << m(i,j);

        if( i != (m.numrows-1) )
            os << std::endl << " ";
    }

    os << " ]";
    os << std::endl;
    return os;
}

}

#include<bigarray.hpp>

#endif // _GURLS_BIGARRAY_H_
