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

#define CHECK_HDF5_ERR(value, errorString) \
    if(value < 0) throw gException(errorString);

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
BigArray<T>::BigArray(): gMat2D<T>(), file_id(-1), dset_id(-1), plist_id(-1), dataFileName("")
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
    close();
    connection.disconnect();
}

template <typename T>
void BigArray<T>::flush()
{
    herr_t status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    CHECK_HDF5_ERR(status, "Error flushing data to file")
}

template <typename T>
void BigArray<T>::close()
{
    herr_t status;

    if(dset_id >= 0)
    {
        status = H5Dclose(dset_id);
        CHECK_HDF5_ERR(status, "Error closing BigArray dataset")
        dset_id = -1;
    }

    if(plist_id >= 0)
    {
        status = H5Pclose(plist_id);
        CHECK_HDF5_ERR(status, "Error closing BigArray plist")
        plist_id = -1;
    }

    if(file_id >= 0)
    {
        status = H5Fclose(file_id);
        CHECK_HDF5_ERR(status, "Error closing BigArray file")
        file_id = -1;
    }
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

    getMatrix(i, 0, 1, this->numcols, ret.getData());

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

    getMatrix(0, i, this->numrows, 1, ret.getData());

    return ret;
}

template <typename T>
gVec<T> BigArray<T>::getColumn(unsigned long i) const
{
    return (*this)(i);
}

template <typename T>
T BigArray<T>::operator() (unsigned long row, unsigned long col) const
{
    if(row >= this->numrows || col >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    T ret;

    std::string errorString("Error reading BigArray value");

    hsize_t dims[2] = {1, 1};
    hid_t memspace = H5Screate_simple(2, dims, NULL);
    CHECK_HDF5_ERR(memspace, errorString)

    hsize_t point[2] = {col, row};
    hid_t filespace = H5Dget_space(dset_id);
    CHECK_HDF5_ERR(filespace, errorString)

    herr_t status;

    status = H5Sselect_elements(filespace, H5S_SELECT_SET, 1, point);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Dread(dset_id, getHdfType<T>(), memspace, filespace, plist_id, &ret);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(memspace);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(filespace);
    CHECK_HDF5_ERR(status, errorString)

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

    getMatrix(startingRow, startingCol, numRows, numCols, result.getData());
}

template <typename T>
void BigArray<T>::getMatrix(unsigned long startingRow, unsigned long startingCol, gMat2D<T>&result) const
{
    if(startingRow >= this->numrows || startingCol >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    if(startingRow+result.rows() > this->numrows || startingCol+result.cols() > this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    getMatrix(startingRow, startingCol, result.rows(), result.cols(), result.getData());
}

template <typename T>
void BigArray<T>::getMatrix(unsigned long startingRow, unsigned long startingCol, unsigned long numRows, unsigned long numCols, T* result) const
{
    std::string errorString("Error reading matrix data");

    hsize_t dims[2] = {numCols, numRows};
    hid_t memspace = H5Screate_simple(2, dims, NULL);
    CHECK_HDF5_ERR(memspace, errorString)

    hsize_t	count[2] = {1, 1};
    hsize_t	stride[2] = {1, 1};
    hsize_t	block[2] = {dims[0], dims[1]};
    hsize_t	offset[2] = {startingCol, startingRow};


    // Select hyperslab in the file.
    hid_t filespace = H5Dget_space(dset_id);
    CHECK_HDF5_ERR(filespace, errorString)

    herr_t status;

    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Dread(dset_id, getHdfType<T>(), memspace, filespace, plist_id, result);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(memspace);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(filespace);
    CHECK_HDF5_ERR(status, errorString)
}

// Setter
template <typename T>
void BigArray<T>::setValue(unsigned long row, unsigned long col, T value)
{
    if(row >= this->numrows || col >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    std::string errorString("Error writing BigArray value");

    hsize_t dims[2] = {1, 1};
    hid_t memspace = H5Screate_simple(2, dims, NULL);
    CHECK_HDF5_ERR(memspace, errorString)

    hsize_t point[2] = {col, row};
    hid_t filespace = H5Dget_space(dset_id);
    CHECK_HDF5_ERR(filespace, errorString)

    herr_t status;

    status = H5Sselect_elements(filespace, H5S_SELECT_SET, 1, point);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Dwrite(dset_id, getHdfType<T>(), memspace, filespace, plist_id, &value);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(memspace);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(filespace);
    CHECK_HDF5_ERR(status, errorString)
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

    std::string errorString("Error writing matrix data");

    hsize_t dims[2] = {M_cols, M_rows};
    hid_t memspace = H5Screate_simple(2, dims, NULL);
    CHECK_HDF5_ERR(memspace, errorString)

    hsize_t	count[2] = {1, 1};
    hsize_t	stride[2] = {1, 1};
    hsize_t	block[2] = {dims[0], dims[1]};
    hsize_t	offset[2] = {startingCol, startingRow};


    // Select hyperslab in the file.
    hid_t filespace = H5Dget_space(dset_id);
    CHECK_HDF5_ERR(filespace, errorString)

    herr_t status;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Dwrite(dset_id, getHdfType<T>(), memspace, filespace, plist_id, M);
    CHECK_HDF5_ERR(status, errorString)

    H5Sclose(memspace);
    CHECK_HDF5_ERR(status, errorString)

    H5Sclose(filespace);
    CHECK_HDF5_ERR(status, errorString)
}

template <typename T>
void BigArray<T>::setColumn(unsigned long col, const gVec<T>& value)
{
    if(col >= this->numcols)
        throw gException(Exception_Index_Out_of_Bound);

    if(value.getSize() != this->numrows)
        throw gException(Exception_Inconsistent_Size);

    setMatrix(0, col, value.getData(), this->numrows, 1);
}

template <typename T>
void BigArray<T>::setRow(unsigned long row, const gVec<T>& value)
{
    if(row >= this->numrows)
        throw gException(Exception_Index_Out_of_Bound);

    if(value.getSize() != this->numcols)
        throw gException(Exception_Inconsistent_Size);

    setMatrix(row, 0, value.getData(), 1, this->numcols);
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

    if(!name.empty())
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

    close();
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

//        flush();
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


    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if(plist_id == -1)
        throw gException(errorString);

    herr_t status;

    status = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    CHECK_HDF5_ERR(status, errorString)

    // Create a new file collectively and release property list identifier.
    file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, plist_id);
    CHECK_HDF5_ERR(file_id, errorString)

    status = H5Pclose(plist_id);
    CHECK_HDF5_ERR(status, errorString)

    dset_id =  H5Dopen(file_id, "mat", H5P_DEFAULT);
    CHECK_HDF5_ERR(dset_id, errorString)

    hid_t filespace = H5Dget_space( dset_id );
    CHECK_HDF5_ERR(filespace, errorString)

    hsize_t dims[2], maxDims[2];
    status = H5Sget_simple_extent_dims(filespace, dims, maxDims);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(filespace);
    CHECK_HDF5_ERR(status, errorString)

    this->numrows = dims[1];
    this->numcols = dims[0];

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(plist_id == -1)
        throw gException(errorString);

    status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    CHECK_HDF5_ERR(status, errorString)

}

template <typename T>
void BigArray<T>::init(std::string& fileName, unsigned long r, unsigned long c)
{
    dataFileName = fileName;

    connection.disconnect();
    connection = releaseSignal().connect(boost::bind(&gurls::BigArray<T>::close, this));

    std::string errorString = "Error creating file " + fileName + ":";

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if(plist_id == -1)
        throw gException(errorString);

    herr_t status;

    status = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    CHECK_HDF5_ERR(status, errorString)

    // Create a new file collectively and release property list identifier.
    file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    CHECK_HDF5_ERR(file_id, errorString)

    status = H5Pclose(plist_id);
    CHECK_HDF5_ERR(status, errorString)


    // Create the dataspace for the dataset.
    hsize_t dims[2];
    dims[0] = static_cast<hsize_t>(c);
    dims[1] = static_cast<hsize_t>(r);
    hid_t filespace = H5Screate_simple(2, dims, NULL);
    CHECK_HDF5_ERR(filespace, errorString)


    hid_t plist_dset_id = H5Pcreate(H5P_DATASET_CREATE);
    if(plist_dset_id == -1)
        throw gException(errorString);

    dset_id = H5Dcreate(file_id, "mat", getHdfType<T>(), filespace, H5P_DEFAULT, plist_dset_id, H5P_DEFAULT);
    CHECK_HDF5_ERR(dset_id, errorString)

    status = H5Pclose(plist_dset_id);
    CHECK_HDF5_ERR(status, errorString)

    status = H5Sclose(filespace);
    CHECK_HDF5_ERR(status, errorString)

    this->numrows = r;
    this->numcols = c;

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(plist_id == -1)
        throw gException(errorString);

    status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    CHECK_HDF5_ERR(status, errorString)

    flush();
}

}

