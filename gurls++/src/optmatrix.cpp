#include "optmatrix.h"

namespace gurls
{

/**
  * OptMatrix empty constructor for float elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<float> >::OptMatrix(): OptMatrixBase () , value (*(new gMat2D<float>()))
{
    this->matType = FLOAT;
}

/**
  * OptMatrix constructor for float elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<float> >::OptMatrix(gMat2D<float>& m): OptMatrixBase(), value(m)
{
    this->matType = FLOAT;
}

/**
  * OptMatrix empty constructor for double elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<double> >::OptMatrix(): OptMatrixBase () , value (*(new gMat2D<double>()))
{
    this->matType = DOUBLE;
}

/**
  * OptMatrix constructor for double elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<double> >::OptMatrix(gMat2D<double>& m): OptMatrixBase(), value(m)
{
    this->matType = DOUBLE;
}

/**
  * OptMatrix empty constructor for unsigned long elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<unsigned long> >::OptMatrix(): OptMatrixBase () , value (*(new gMat2D<unsigned long>()))
{
    this->matType = ULONG;
}

/**
  * OptMatrix constructor for unsigned long elements
  */
template <>
GURLS_EXPORT OptMatrix <gMat2D<unsigned long> >::OptMatrix(gMat2D<unsigned long>& m): OptMatrixBase(), value(m)
{
    this->matType = ULONG;
}

}
