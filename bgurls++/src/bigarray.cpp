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
    return netCDF::NcUint();
}

}
