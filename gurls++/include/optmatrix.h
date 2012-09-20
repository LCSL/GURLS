/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, IIT@MIT Lab
 * All rights reserved.
 *
 * authors:  M. Santoro
 * email:   msantoro@mit.edu
 * website: http://cbcl.mit.edu/IIT@MIT/IIT@MIT.html
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name(s) of the copyright holders nor the names
 *       of its contributors or of the Massacusetts Institute of
 *       Technology or of the Italian Institute of Technology may be
 *       used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef _GURLS_OPTMATRIX_H_
#define _GURLS_OPTMATRIX_H_

#include <options.h>

namespace gurls
{

/**
  * \ingroup Settings
  * \brief OptMatrixBase is the base class for all options containing matrices.
  */
class GURLS_EXPORT OptMatrixBase: public GurlsOption
{
public:

    /**
      * Empty constructor
      */
    OptMatrixBase(): GurlsOption(MatrixOption){}

    /**
      * \enum MatrixType
      * Enumeration containing all supported element types
      */
    enum MatrixType{FLOAT, DOUBLE, ULONG};

    /**
      * Returns the element type for the matrix
      */
    MatrixType getMatrixType() const
    {
        return matType;
    }

protected:
    MatrixType matType; ///< Stores the type of the elements inside the matrix
};

/**
  * \ingroup Settings
  * \brief OptMatrix is an option containing a matrix.
  * \tparam MatrixType Type of the matrix contained into the option
  */
template <typename Matrix>
class OptMatrix: public OptMatrixBase
{
private:
    Matrix& value;  ///< Option value

public:
    typedef Matrix ValueType;

    /**
      * Empty constructor
      */
    OptMatrix(): OptMatrixBase () , value (*(new Matrix(2,2)))
    {
        throw gException(Exception_Unsupported_MatrixType);
    }

    /**
      * Constructor from an existing matrix
      */
    OptMatrix(Matrix& m): OptMatrixBase(), value(m)
    {
        throw gException(Exception_Unsupported_MatrixType);
    }

    /**
      * Copies the option values from an existing \ref OptMatrix
      */
    OptMatrix<Matrix>& operator=(const OptMatrix<Matrix>& other);

    /**
      * Destructor
      */
    ~OptMatrix()
    {
        delete &value;
    }

    /**
      * Copies the matrix from an existing matrix
      */
    OptMatrix& operator=(const Matrix& other)
    {
        this->type = MatrixOption;
        this->value = other;

        return *this;
    }

    /**
      * Copies the matrix from an existing matrix
      */
    void setValue(const Matrix& newvalue)
    {
        value = newvalue;
    }

    /**
      * Returns the matrix
      */
    Matrix& getValue()
    {
        return value;
    }

    /**
      * Returns the matrix
      */
    const Matrix& getValue() const
    {
        return value;
    }

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const
    {
        return (id == MatrixOption);
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptMatrix
      */
    static OptMatrix* dynacast(GurlsOption* opt)
    {
        if (opt->isA(MatrixOption) )
            return static_cast<OptMatrix*>(opt);

        throw gException(gurls::Exception_Illegal_Dynamic_Cast);
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref OptMatrix
      */
    static const OptMatrix* dynacast(const GurlsOption* opt)
    {
        if (opt->isA(MatrixOption) )
            return static_cast<const OptMatrix*>(opt);

        throw gException(gurls::Exception_Illegal_Dynamic_Cast);
    }

    /**
      * Writes the option to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

};

/**
  * OptMatrix empty constructor for float elements
  */
template <>
OptMatrix <gMat2D<float> >::OptMatrix();

/**
  * OptMatrix constructor for float elements
  */
template <>
OptMatrix <gMat2D<float> >::OptMatrix(gMat2D<float>& m);

/**
  * OptMatrix empty constructor for double elements
  */
template <>
OptMatrix <gMat2D<double> >::OptMatrix();

/**
  * OptMatrix constructor for double elements
  */
template <>
OptMatrix <gMat2D<double> >::OptMatrix(gMat2D<double>& m);

/**
  * OptMatrix empty constructor for unsigned long elements
  */
template <>
OptMatrix <gMat2D<unsigned long> >::OptMatrix();

/**
  * OptMatrix constructor for unsigned long elements
  */
template <>
OptMatrix <gMat2D<unsigned long> >::OptMatrix(gMat2D<unsigned long>& m);


/**
  * Writes an OptMatrix to a stream
  */
template <typename T>
std::ostream& OptMatrix<T>::operator << (std::ostream& os)
{
    return os << std::endl << this->getValue();
}

}

#endif //_GURLS_OPTMATRIX_H_
