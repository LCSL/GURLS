#include "optmatrix.h"

namespace gurls
{

template<>
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<gMat2D<float> >()
{
    return OptMatrixBase::FLOAT;
}

template<>
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<gMat2D<double> >()
{
    return OptMatrixBase::DOUBLE;
}

template<>
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<gMat2D<unsigned long> >()
{
    return OptMatrixBase::ULONG;
}

#ifdef _BGURLS

template<>
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<BigArray<float> >()
{
    return OptMatrixBase::FLOAT;
}

template<>
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<BigArray<double> >()
{
    return OptMatrixBase::DOUBLE;
}

template<>
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<BigArray<unsigned long> >()
{
    return OptMatrixBase::ULONG;
}

template <>
GURLS_EXPORT bool OptMatrix <BigArray<float> >::hasBigArray() const
{
    return true;
}

template <>
GURLS_EXPORT bool OptMatrix <BigArray<double> >::hasBigArray() const
{
    return true;
}

template <>
GURLS_EXPORT bool OptMatrix <BigArray<unsigned long> >::hasBigArray() const
{
    return true;
}

#endif

}
