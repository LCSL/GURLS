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

}
