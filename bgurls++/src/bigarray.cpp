#include "bgurls++/bigarray.h"

namespace gurls
{

template<>
GURLS_EXPORT hid_t getHdfType<float>()
{
    return H5T_NATIVE_FLOAT;
}

template<>
GURLS_EXPORT  hid_t getHdfType<double>()
{
    return H5T_NATIVE_DOUBLE;
}

template<>
GURLS_EXPORT hid_t getHdfType<unsigned long>()
{
    return H5T_NATIVE_ULONG;
}

template<>
GURLS_EXPORT hid_t getHdfType<unsigned int>()
{
    return H5T_NATIVE_UINT;
}

}
