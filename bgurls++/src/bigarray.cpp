#include "bigarray.h"

namespace gurls
{

template<>
GURLS_EXPORT nc_type getNcType<float>()
{
    return NC_FLOAT;
}

template<>
GURLS_EXPORT  nc_type getNcType<double>()
{
    return NC_DOUBLE;
}

template<>
GURLS_EXPORT nc_type getNcType<unsigned long>()
{
    if(sizeof(unsigned long) == sizeof(unsigned int))
        return NC_UINT;

    return NC_UINT64;
}

template<>
GURLS_EXPORT nc_type getNcType<unsigned int>()
{
    return NC_UINT;
}


template<>
GURLS_EXPORT int nc_put_vara<float>(int ncid, int varid,  const size_t *startp,
                                    const size_t *countp, const float *op)
{
    return nc_put_vara_float(ncid, varid, startp, countp, op);
}

template<>
GURLS_EXPORT int nc_put_vara<double>(int ncid, int varid,  const size_t *startp,
                                     const size_t *countp, const double *op)
{
    return nc_put_vara_double(ncid, varid, startp, countp, op);
}

template<>
GURLS_EXPORT int nc_put_vara<unsigned long>(int ncid, int varid,  const size_t *startp,
                                            const size_t *countp, const unsigned long *op)
{
    if(sizeof(unsigned long) == sizeof(unsigned int))
        return nc_put_vara_uint(ncid, varid, startp, countp, (const unsigned int*) op);

    return nc_put_vara_ulonglong(ncid, varid, startp, countp, (const unsigned long long*) op);
}

template<>
GURLS_EXPORT int nc_put_vara<unsigned int>(int ncid, int varid,  const size_t *startp,
                                           const size_t *countp, const unsigned int *op)
{
    return nc_put_vara_uint(ncid, varid, startp, countp, op);
}



template<>
GURLS_EXPORT int nc_get_vara<float>(int ncid, int varid,  const size_t *startp,
                                    const size_t *countp, float *op)
{
    return nc_get_vara_float(ncid, varid, startp, countp, op);
}

template<>
GURLS_EXPORT int nc_get_vara<double>(int ncid, int varid,  const size_t *startp,
                                     const size_t *countp, double *op)
{
    return nc_get_vara_double(ncid, varid, startp, countp, op);
}

template<>
GURLS_EXPORT int nc_get_vara<unsigned long>(int ncid, int varid,  const size_t *startp,
                                            const size_t *countp, unsigned long *op)
{
    if(sizeof(unsigned long) == sizeof(unsigned int))
        return nc_get_vara_uint(ncid, varid, startp, countp, (unsigned int*) op);

    return nc_get_vara_ulonglong(ncid, varid, startp, countp, (unsigned long long*) op);
}

template<>
GURLS_EXPORT int nc_get_vara<unsigned int>(int ncid, int varid,  const size_t *startp,
                                           const size_t *countp, unsigned int *op)
{
    return nc_get_vara_uint(ncid, varid, startp, countp, op);
}


template<>
GURLS_EXPORT int nc_get_var1<float>(int ncid, int varid,  const size_t *indexp, float *op)
{
    return nc_get_var1_float(ncid, varid, indexp, op);
}

template<>
GURLS_EXPORT int nc_get_var1<double>(int ncid, int varid,  const size_t *indexp, double *op)
{
    return nc_get_var1_double(ncid, varid, indexp, op);
}

template<>
GURLS_EXPORT int nc_get_var1<unsigned long>(int ncid, int varid,  const size_t *indexp, unsigned long *op)
{
    if(sizeof(unsigned long) == sizeof(unsigned int))
        return nc_get_var1_uint(ncid, varid, indexp, (unsigned int*) op);

    return nc_get_var1_ulonglong(ncid, varid, indexp, (unsigned long long*) op);
}

template<>
GURLS_EXPORT int nc_get_var1<unsigned int>(int ncid, int varid,  const size_t *indexp, unsigned int *op)
{
    return nc_get_var1_uint(ncid, varid, indexp, op);
}


template<>
GURLS_EXPORT int nc_put_var1<float>(int ncid, int varid,  const size_t *indexp, const float *op)
{
    return nc_put_var1_float(ncid, varid, indexp, op);
}

template<>
GURLS_EXPORT int nc_put_var1<double>(int ncid, int varid,  const size_t *indexp, const double *op)
{
    return nc_put_var1_double(ncid, varid, indexp, op);
}

template<>
GURLS_EXPORT int nc_put_var1<unsigned long>(int ncid, int varid,  const size_t *indexp, const unsigned long *op)
{
    if(sizeof(unsigned long) == sizeof(unsigned int))
        return nc_put_var1_uint(ncid, varid, indexp, (unsigned int*) op);

    return nc_put_var1_ulonglong(ncid, varid, indexp, (unsigned long long*) op);
}

template<>
GURLS_EXPORT int nc_put_var1<unsigned int>(int ncid, int varid,  const size_t *indexp, const unsigned int *op)
{
    return nc_put_var1_uint(ncid, varid, indexp, op);
}

}
