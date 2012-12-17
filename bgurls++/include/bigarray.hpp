#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/base_object.hpp>


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

}
