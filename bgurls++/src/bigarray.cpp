#include "bigarray.h"

namespace gurls
{

template<>
GURLS_EXPORT netCDF::NcType getNcType<float>()
{
    return netCDF::NcFloat();
}

template<>
GURLS_EXPORT  netCDF::NcType getNcType<double>()
{
    return netCDF::NcDouble();
}

template<>
GURLS_EXPORT netCDF::NcType getNcType<unsigned long>()
{
    if(sizeof(unsigned long) == sizeof(unsigned int))
        return netCDF::NcUint();

    return netCDF::NcType(NC_UINT64);
}

template<>
GURLS_EXPORT netCDF::NcType getNcType<unsigned int>()
{
    return netCDF::NcUint();
}

}
