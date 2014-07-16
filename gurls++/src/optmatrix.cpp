#include "gurls++/optmatrix.h"

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
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<const gMat2D<float> >()
{
    return OptMatrixBase::FLOAT;
}

template<>
GURLS_EXPORT OptMatrixBase::MatrixType getMatrixCellType<const gMat2D<double> >()
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

template<>
GURLS_EXPORT bool containsBigArray<BigArray<float> >()
{
    return true;
}


template<>
GURLS_EXPORT bool containsBigArray<BigArray<double> >()
{
    return true;
}


template<>
GURLS_EXPORT bool containsBigArray<BigArray<unsigned long> >()
{
    return true;
}


#endif

}
