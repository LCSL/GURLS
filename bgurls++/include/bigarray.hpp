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

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/bind.hpp>

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
BigArray<T>::BigArray(): gMat2D<T>(), dataFile(0), dataFileName("")
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

    MPI_Barrier(MPI_COMM_WORLD);
}

template <typename T>
BigArray<T>::BigArray(std::string fileName, const gMat2D<T>& mat)
{
    init(fileName, mat.rows(), mat.cols());

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid ==  0)
        setMatrix(0, 0, mat);

    MPI_Barrier(MPI_COMM_WORLD);
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
    nc_close(dataFile);
    connection.disconnect();
}

template <typename T>
void BigArray<T>::flush()
{
//    nc_close(dataFile);
//    loadNC(dataFileName);
}

template <typename T>
void BigArray<T>::close()
{

    nc_close(dataFile);
    dataFile = 0;
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

    size_t start[2], count[2];

    start[0] = 0;
    start[1] = i;

    count[0] = this->numcols;
    count[1] = 1;

    int retval;
    if ((retval = nc_get_vara<T>(dataFile, matrix, start, count, ret.getData())))
       throw gException("Error reading row");

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

    size_t start[2], count[2];

    start[0] = i;
    start[1] = 0;

    count[0] = 1;
    count[1] = this->numrows;

    int retval;
    if ((retval = nc_get_vara<T>(dataFile, matrix, start, count, ret.getData())))
       throw gException("Error reading column");

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

    size_t start[2];

    start[0] = col;
    start[1] = row;

    int retval;
    if ((retval = nc_get_var1<T>(dataFile, matrix, start, &ret)))
       throw gException("Error reading value");

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

    size_t start[2], count[2];

    start[0] = startingCol;
    start[1] = startingRow;

    count[0] = numCols;
    count[1] = numRows;

    int retval;
    if ((retval = nc_get_vara<T>(dataFile, matrix, start, count, result.getData())))
       throw gException("Error reading matrix");

}

template <typename T>
void BigArray<T>::getMatrix(unsigned long startingRow, unsigned long startingCol, gMat2D<T>&result) const
{
    if(startingRow >= this->numrows || startingCol >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    if(startingRow+result.rows() > this->numrows || startingCol+result.cols() > this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    size_t start[2], count[2];

    start[0] = startingCol;
    start[1] = startingRow;

    count[0] = result.cols();
    count[1] = result.rows();

    int retval;
    if ((retval = nc_get_vara<T>(dataFile, matrix, start, count, result.getData())))
       throw gException("Error reading matrix");
}



// Setter
template <typename T>
void BigArray<T>::setValue(unsigned long row, unsigned long col, T value)
{
    if(row >= this->numrows || col >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    size_t start[2];

    start[0] = col;
    start[1] = row;

    int retval;
    if ((retval = nc_put_var1<T>(dataFile, matrix, start, &value)))
        throw gException("Error setting value");
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

    size_t start[2], count[2];

    start[0] = startingCol;
    start[1] = startingRow;

    count[0] = M_cols;
    count[1] = M_rows;

    int retval;
    if ((retval = nc_put_vara<T>(dataFile, matrix, start, count, M)))
        throw gException("Error setting matrix");
}

template <typename T>
void BigArray<T>::setColumn(unsigned long col, const gVec<T>& value)
{
    if(col >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    if(value.getSize() != this->numrows)
        throw gException(Exception_Inconsistent_Size);

    size_t start[2], count[2];

    start[0] = col;
    start[1] = 0;

    count[0] = 1;
    count[1] = this->numrows;

    int retval;
    if ((retval = nc_put_vara<T>(dataFile, matrix, start, count, value.getData())))
        throw gException("Error setting column");
}

template <typename T>
void BigArray<T>::setRow(unsigned long row, const gVec<T>& value)
{
    if(row >= this->numrows)
        throw gException(Exception_Index_Out_of_Bound);

    if(value.getSize() != this->numcols)
        throw gException(Exception_Inconsistent_Size);

    size_t start[2], count[2];

    start[0] = 0;
    start[1] = row;

    count[0] = this->numcols;
    count[1] = 1;

    int retval;
    if ((retval = nc_put_vara<T>(dataFile, matrix, start, count, value.getData())))
        throw gException("Error setting row");

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
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    unsigned long rows = 0;
    unsigned long cols = 0;

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(" ;|,");
    std::string line;

    if(myid == 0)
    {
        std::ifstream in(fileName.c_str());

        if(!in.is_open())
            throw gurls::gException("Cannot open file " + fileName);

        while(std::getline(in, line))
        {
            if(!line.empty())
                ++rows;

            if(rows == 1)
            {
                tokenizer tokens(line, sep);
                for (tokenizer::iterator t_it = tokens.begin(); t_it != tokens.end(); ++t_it)
                    ++cols;
            }
        }

        in.close();
    }

    MPI_Bcast(&rows, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    nc_close(dataFile);
    init(dataFileName, rows, cols);

    if(rows == 0 || cols == 0)
        return;

    if(myid == 0)
    {
        std::ifstream in(fileName.c_str());

        if(!in.is_open())
            throw gurls::gException("Cannot open file " + fileName);

        gVec<T> row_vector(cols);
        T* rv_it = NULL;
        T* rv_end = row_vector.getData()+cols;
        tokenizer::iterator t_it;

        unsigned long row = 0;
        while (std::getline(in, line))
        {
            if(!line.empty())
            {
                tokenizer tokens(line, sep);
                for (rv_it = row_vector.getData(), t_it = tokens.begin(); rv_it != rv_end; ++t_it, ++rv_it)
                        *rv_it = boost::lexical_cast<T>(*t_it);

                if(t_it != tokens.end())
                    throw gException(Exception_Inconsistent_Size);

                setRow(row, row_vector);

                ++row;
            }
        }
        in.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}


template <typename T>
void BigArray<T>::loadNC(const std::string &fileName)
{
    dataFileName = fileName;

    connection.disconnect();
    connection = releaseSignal().connect(boost::bind(&gurls::BigArray<T>::close, this));

    std::string errorString = "Error opening file " + fileName + ":";
    int retval;
    char buf[NC_MAX_NAME+1];

    if((retval = nc_open_par(fileName.c_str(), NC_MPIIO|NC_WRITE, MPI_COMM_WORLD, MPI_INFO_NULL, &dataFile)))
        throw gException(errorString);


    int dims[2];
    if((retval = nc_inq_dimid(dataFile, "rows", &dims[0])))
        throw gException(errorString);
    if((retval = nc_inq_dimid(dataFile, "cols", &dims[1])))
        throw gException(errorString);


    if((retval = nc_inq_dim(dataFile, dims[0], buf, &(this->numcols))))
        throw gException(errorString);
    if((retval = nc_inq_dim(dataFile, dims[1], buf, &(this->numrows))))
        throw gException(errorString);


    if((retval = nc_inq_varid(dataFile, "mat", &matrix)))
        throw gException(errorString);

}

template <typename T>
void BigArray<T>::init(std::string& fileName, unsigned long r, unsigned long c)
{
    dataFileName = fileName;

    connection.disconnect();
    connection = releaseSignal().connect(boost::bind(&gurls::BigArray<T>::close, this));

    std::string errorString = "Error opening file " + fileName + ":";
    int retval;

    if((retval = nc_create_par(fileName.c_str(), NC_MPIIO|NC_NETCDF4|NC_CLOBBER,  MPI_COMM_WORLD, MPI_INFO_NULL, &dataFile)))
        throw gException(errorString);

    int dims[2];

    if ((retval = nc_def_dim(dataFile, "rows", c, &dims[0])))
       throw gException(errorString);
    if ((retval = nc_def_dim(dataFile, "cols", r, &dims[1])))
       throw gException(errorString);

    this->numrows = r;
    this->numcols = c;

    if ((retval = nc_def_var(dataFile, "mat", getNcType<T>(), 2, dims, &this->matrix)))
        throw gException(errorString);

    if(nc_enddef(dataFile))
        throw gException(errorString);

}

}
