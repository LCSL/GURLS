#include "optarray.h"
#include "serialization.h"

#include <fstream>

using namespace std;

namespace gurls{

OptArray::OptArray():GurlsOption(OptArrayOption)
{
    value = new ValueType();
}

OptArray::OptArray(const OptArray& other):GurlsOption(OptArrayOption)
{
    value = new ValueType(*(other.value));
}

OptArray::~OptArray()
{
    clear();

    delete value;
}

const OptArray::ValueType& OptArray::getValue() const
{
    return *value;
}

void OptArray::push_back(GurlsOption* opt)
{
    value->push_back(opt);
}

void OptArray::erase(unsigned long index, bool deleteMembers)
{
    GurlsOption* opt = (*this)[index];

    value->erase(value->begin()+index);

    if(deleteMembers)
        delete opt;
}

void OptArray::clear()
{
    for (ValueType::iterator it = value->begin(); it != value->end(); ++it)
        delete *it;

    value->clear();
}

void OptArray::reserve(unsigned long size)
{
    value->reserve(size);
}

bool OptArray::isA(OptTypes id) const
{
    return (id == OptArrayOption);
}

unsigned long OptArray::size() const
{
    return value->size();
}

OptArray* OptArray::dynacast(GurlsOption* opt)
{
    if (opt->isA(OptArrayOption))
        return static_cast<OptArray*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

const OptArray* OptArray::dynacast(const GurlsOption* opt)
{
    if (opt->isA(OptArrayOption))
        return static_cast<const OptArray*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

GurlsOption* OptArray::operator[] (unsigned long i)
{
    if ( i >= value->size() )
        throw gException(gurls::Exception_Index_Out_of_Bound);

    return (*value)[i];
}

/**
  * Writes an OptArray to a stream
  */
GURLS_EXPORT std::ostream& operator<<(std::ostream& os, OptArray& opt)
{
    os << endl << "~~~~~~~ OptArray: " << endl;

    const unsigned long size = opt.size();

    for(unsigned long i=0; i< size; ++i)
        os << "[ " << i << " ] = " << (*(opt[i])) << endl;

    os << "~~~~~~~";

    return os;
}

std::ostream& OptArray::operator<<(std::ostream& os)
{
    return os << *this;
}

void OptArray::save(const std::string& fileName) const
{
#ifndef USE_BINARY_ARCHIVES
    std::ofstream outstream(fileName.c_str());
#else
    std::ofstream outstream(fileName.c_str(), ios_base::binary);
#endif

    if(!outstream.is_open())
        throw gException("Could not open file " + fileName);

    oarchive outar(outstream);
    outar << *this;

    outstream.close();
}


void OptArray::load(const std::string& fileName)
{
#ifndef USE_BINARY_ARCHIVES
    std::ifstream instream(fileName.c_str());
#else
    std::ifstream instream(fileName.c_str(), ios_base::binary);
#endif

    if(!instream.is_open())
        throw gException("Could not open file " + fileName);

    try
    {
        iarchive inar(instream);
        inar >> *this;
    }
    catch(boost::archive::archive_exception&)
    {
        instream.close();
        throw gException("Invalid file format for " + fileName);
    }

    instream.close();
}

}
