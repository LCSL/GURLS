#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/base_object.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <exceptions.h>

#ifdef  USE_BINARY_ARCHIVES
typedef boost::archive::binary_iarchive iarchive;
typedef boost::archive::binary_oarchive oarchive;
#else
typedef boost::archive::text_iarchive iarchive;
typedef boost::archive::text_oarchive oarchive;
#endif

namespace gurls
{

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

template <typename T>
BigArray<T>::BigArray(): gMat2D<T>(), dataFile(NULL), dataFileName("")
{
}

template <typename T>
BigArray<T>::BigArray(std::string fileName)
{
    loadNC(fileName);
}

template <typename T>
BigArray<T>::BigArray(std::string fileName, unsigned long r, unsigned long c)
{
    init(fileName, r, c);

    if(r>0 && c>0)
        setValue(0, 0, 0);

    flush();
}

template <typename T>
BigArray<T>::BigArray(std::string fileName, const gMat2D<T>& mat)
{
    init(fileName, mat.rows(), mat.cols());
    setMatrix(0, 0, mat);
    flush();
}

template <typename T>
BigArray<T>::BigArray(const BigArray<T>& other)
{
    loadNC(other.dataFileName);
}

//template <typename T>
//BigArray<T>::BigArray(const BigArray<T>& other, std::string fileName){}

//template <typename T>
//BigArray<T>::BigArray<T>& operator= (BigArray<T> other)
//{
//    if(!dataFileName.empty())
//        delete dataFile;
//    loadNC(other.dataFileName);
//    return *this;
//}

template <typename T>
BigArray<T>::~BigArray()
{
    delete dataFile;
}

template <typename T>
void BigArray<T>::flush()
{
    delete dataFile;
    loadNC(dataFileName);
}

template <typename T>
unsigned long BigArray<T>::rows() const
{
    return this->numrows;
}

template <typename T>
unsigned long BigArray<T>::cols() const
{
    return this->numcols;
}


template <typename T>
std::string BigArray<T>::what() const
{
    std::stringstream v;
    v << "BigArray: (" << this->numrows;
    v << " x " << this->numcols;
    v << ") matrix of type " << typeid(T).name();
    return v.str();
}

// Getter
template <typename T>
gVec<T> BigArray<T>::operator[](unsigned long i) const
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

template <typename T>
gVec<T> BigArray<T>::getRow(unsigned long i)
{
    return (*this)[i];
}

template <typename T>
gVec<T> BigArray<T>::operator() (unsigned long i) const
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

template <typename T>
gVec<T> BigArray<T>::getColumnn(unsigned long i) const
{
    return (*this)(i);
}

template <typename T>
T BigArray<T>::operator() (unsigned long row, unsigned long col) const
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

template <typename T>
T BigArray<T>::getValue(unsigned long row, unsigned long col) const
{
    return (*this)(row, col);
}

template <typename T>
void BigArray<T>::getMatrix(unsigned long startingRow, unsigned long startingCol, unsigned long numRows, unsigned long numCols, gMat2D<T>&result) const
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

template <typename T>
void BigArray<T>::getMatrix(unsigned long startingRow, unsigned long startingCol, gMat2D<T>&result) const
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
template <typename T>
void BigArray<T>::setValue(unsigned long row, unsigned long col, T value)
{
    if(row >= this->numrows || col >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    std::vector<size_t> start(2), count(2, 1);

    start[0] = col;
    start[1] = row;

    matrix.putVar(start, count, &value);
}

template <typename T>
void BigArray<T>::setMatrix(unsigned long startingRow, unsigned long startingCol, const gMat2D<T>&value)
{
    setMatrix(startingRow, startingCol, value.getData(), value.rows(), value.cols());
}

template <typename T>
void BigArray<T>::setMatrix(unsigned long startingRow, unsigned long startingCol, const T* M, const unsigned long M_rows, const unsigned long M_cols)
{
    if(startingRow >= this->numrows || startingCol >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    if(startingRow+M_rows > this->numrows || startingCol+M_cols > this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    std::vector<size_t> start(2), count(2);

    start[0] = startingCol;
    start[1] = startingRow;

    count[0] = M_cols;
    count[1] = M_rows;

    matrix.putVar(start, count, M);
}

template <typename T>
void BigArray<T>::setColumn(unsigned long col, const gVec<T>& value)
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

template <typename T>
void BigArray<T>::setRow(unsigned long row, const gVec<T>& value)
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


template <typename T>
template<class Archive>
void BigArray<T>::save(Archive & ar, const unsigned int /* file_version */) const
{
//    ar & this->dataFileName & this->numrows & this->numcols;

//    T value;
//    for(unsigned long i=0; i< this->numrows; ++i)
//    {
//        for(unsigned long j=0; j< this->numcols; ++j)
//        {
//            value = (*this)(i,j);
//            ar & value;
//        }
//    }

    ar & this->dataFileName;

}

template <typename T>
template<class Archive>
void BigArray<T>::load(Archive & ar, const unsigned int /* file_version */)
{
//    std::string name;
//    unsigned long rows, cols;

//    ar & name & rows & cols;

//    name = name + "_deserialized";
//    init(name, rows, cols);

//    T value;
//    for(unsigned long i=0; i< this->numrows; ++i)
//    {
//        for(unsigned long j=0; j< this->numcols; ++j)
//        {
//            ar & value;
//            setValue(i, j, value);
//        }
//    }

    std::string name;
    ar & name;

    loadNC(name);
}

template <typename T>
void BigArray<T>::readCSV(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if(!in.is_open())
        throw gurls::gException("Cannot open file " + fileName);

    unsigned long rows = 0;
    unsigned long cols = 0;

    std::string line;

    while(std::getline(in, line))
    {
        if(!line.empty())
            ++rows;

        if(rows == 1)
        {
            std::istringstream ss(line);
            std::vector<T> tf;
            std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(), std::back_inserter(tf));

            cols = tf.size();
        }
    }

    in.clear();
    in.seekg (0, std::ios_base::beg);

    delete dataFile;
    init(dataFileName, rows, cols);

    if(rows == 0 || cols == 0)
        return;

    unsigned long row = 0;
    while (std::getline(in, line))
    {
        if(!line.empty())
        {
            std::istringstream ss(line);
            std::vector<T> tf;
            std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(), std::back_inserter(tf));

            if(tf.size() != cols)
                throw gException(Exception_Inconsistent_Size);

            for(unsigned long i=0; i<cols; ++i)
                setValue(row, i, tf[i]);

            ++row;
        }
    }

    in.close();
}


template <typename T>
void BigArray<T>::loadNC(const std::string &fileName)
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

template <typename T>
void BigArray<T>::init(std::string& fileName, unsigned long r, unsigned long c)
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

}
