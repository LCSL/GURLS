/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011-1013, IIT@MIT Lab
  * All rights reserved.
  *
  * author:  M. Santoro
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

#ifndef _GURLS_OPTSERIALIZATION_H_
#define _GURLS_OPTSERIALIZATION_H_

#include <string>

#include "optarray.h"
#include "optlist.h"
#include "optfunction.h"
#include "optmatrix.h"

#include <boost/serialization/split_free.hpp>

namespace boost {
namespace serialization {

/**
  * Serializes an OptArray to a generic archive
  */
template<class Archive>
void save(Archive & ar, const gurls::OptArray & opt, const unsigned int version)
{
    gurls::OptTypes type;

    type = opt.getType();
    ar & type;

    unsigned long size = opt.size();
    ar & size;

    const gurls::OptArray::ValueType& value = opt.getValue();

    for (gurls::OptArray::ValueType::const_iterator it = value.begin(); it != value.end(); ++it)
    {
        gurls::GurlsOption* opt0 = *it;
        type = opt0->getType();

        ar & type;

        switch(type)
        {
        case gurls::StringOption:
            ar & (*static_cast<gurls::OptString*>(opt0));
            break;
        case gurls::StringListOption:
            ar & (*static_cast<gurls::OptStringList*>(opt0));
            break;
        case gurls::NumberOption:
            ar & (*static_cast<gurls::OptNumber*>(opt0));
            break;
        case gurls::NumberListOption:
            ar & (*static_cast<gurls::OptNumberList*>(opt0));
            break;
        case gurls::FunctionOption:
            ar & (*static_cast<gurls::OptFunction*>(opt0));
            break;
        case gurls::MatrixOption:
        {
            gurls::OptMatrixBase* optM = static_cast<gurls::OptMatrixBase*>(opt0);
            gurls::OptMatrixBase::MatrixType matType = optM->getMatrixType();

            ar & matType;

#ifdef _BGURLS

            bool bigArray = optM->hasBigArray();
            ar & bigArray;

            if(bigArray)
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<float> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<double> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<unsigned long> >*>(optM));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }
            else
#endif
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<float> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<double> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<unsigned long> >*>(optM));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }
            break;

        }
        case gurls::TaskSequenceOption:
            ar & (*static_cast<gurls::OptTaskSequence*>(opt0));
            break;
        case gurls::ProcessOption:
            ar & (*static_cast<gurls::OptProcess*>(opt0));
            break;
        case gurls::OptListOption:
            ar & (*static_cast<gurls::GurlsOptionsList*>(opt0));
            break;
        case gurls::OptArrayOption:
            ar & (*static_cast<gurls::OptArray*>(opt0));
            break;
        default:
            throw gurls::gException(gurls::Exception_Unknown_Option);
        }
    }
}

/**
  * Deserializes an OptArray from a generic archive
  */
template<class Archive>
void load(Archive &ar, gurls::OptArray &t, const unsigned int)
{
    gurls::OptTypes type;

    type = t.getType();
    ar & type;

    unsigned long n = 0;
    ar & n;

    t.clear();
    t.reserve(n);

    for (unsigned long i = 0; i<n; ++i)
    {
        ar & type;

        switch(type)
        {

        case gurls::StringOption:
        {
            gurls::OptString* opt = new gurls::OptString();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::StringListOption:
        {
            gurls::OptStringList* opt = new gurls::OptStringList();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::NumberOption:
        {
            gurls::OptNumber* opt = new gurls::OptNumber();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::NumberListOption:
        {
            gurls::OptNumberList* opt = new gurls::OptNumberList();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::FunctionOption:
        {
            gurls::OptFunction* opt = new gurls::OptFunction();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::MatrixOption:
        {

            gurls::OptMatrixBase::MatrixType matType;
            ar & matType;

            gurls::OptMatrixBase* opt;

#ifdef _BGURLS

            bool bigArray;
            ar & bigArray;

            if(bigArray)
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        opt = new gurls::OptMatrix<gurls::BigArray<float> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<float> >*>(opt));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        opt = new gurls::OptMatrix<gurls::BigArray<double> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<double> >*>(opt));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        opt = new gurls::OptMatrix<gurls::BigArray<unsigned long> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<unsigned long> >*>(opt));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }
            else
#endif
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        opt = new gurls::OptMatrix<gurls::gMat2D<float> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<float> >*>(opt));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        opt = new gurls::OptMatrix<gurls::gMat2D<double> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<double> >*>(opt));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        opt = new gurls::OptMatrix<gurls::gMat2D<unsigned long> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<unsigned long> >*>(opt));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }

            t.push_back(opt);
            break;

        }
        case gurls::TaskSequenceOption:
        {
            gurls::OptTaskSequence* opt = new gurls::OptTaskSequence();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::ProcessOption:
        {
            gurls::OptProcess* opt = new gurls::OptProcess();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::OptListOption:
        {
            gurls::GurlsOptionsList* opt = new gurls::GurlsOptionsList("", false);
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        case gurls::OptArrayOption:
        {
            gurls::OptArray* opt = new gurls::OptArray();
            ar & (*opt);
            t.push_back(opt);
            break;
        }
        default:
            throw gurls::gException(gurls::Exception_Unknown_Option);

        }
    }
}



/**
  * Serializes a GurlsOptionsList to a generic archive
  */
template<class Archive>
void save(Archive & ar, const gurls::GurlsOptionsList& opt, const unsigned int /* file_version */)
{
    gurls::OptTypes type;

    type = opt.getType();
    ar & type;

    std::string name = opt.getName();
    ar & name;

    int n = opt.size();
    ar & n;

    typedef gurls::GurlsOptionsList::ValueType ValueType;
    const ValueType& value = opt.getValue();

    for (ValueType::const_iterator it = value.begin(); it != value.end(); ++it)
    {
        name = it->first;
        ar & name;

        gurls::GurlsOption* opt0 = it->second;
        type = opt0->getType();
        ar & type;


        switch(type)
        {
        case gurls::StringOption:
            ar & (*static_cast<gurls::OptString*>(opt0));
            break;
        case gurls::StringListOption:
            ar & (*static_cast<gurls::OptStringList*>(opt0));
            break;
        case gurls::NumberOption:
            ar & (*static_cast<gurls::OptNumber*>(opt0));
            break;
        case gurls::NumberListOption:
            ar & (*static_cast<gurls::OptNumberList*>(opt0));
            break;
        case gurls::FunctionOption:
            ar & (*static_cast<gurls::OptFunction*>(opt0));
            break;
        case gurls::MatrixOption:
        {
            gurls::OptMatrixBase* optM = static_cast<gurls::OptMatrixBase*>(opt0);
            gurls::OptMatrixBase::MatrixType matType = optM->getMatrixType();

            ar & matType;

#ifdef _BGURLS

            bool bigArray = optM->hasBigArray();
            ar & bigArray;

            if(bigArray)
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<float> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<double> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<unsigned long> >*>(optM));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }
            else
#endif
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<float> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<double> >*>(optM));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<unsigned long> >*>(optM));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }

            break;

        }
        case gurls::TaskSequenceOption:
            ar & (*static_cast<gurls::OptTaskSequence*>(opt0));
            break;
        case gurls::ProcessOption:
            ar & (*static_cast<gurls::OptProcess*>(opt0));
            break;
        case gurls::OptListOption:
            ar & (*static_cast<gurls::GurlsOptionsList*>(opt0));
            break;
        case gurls::OptArrayOption:
            ar & (*static_cast<gurls::OptArray*>(opt0));
            break;
        default:
            throw gurls::gException(gurls::Exception_Unknown_Option);
        }
    }
}

/**
  * Deserializes a GurlsOptionsList from a generic archive
  */
template<class Archive>
void load(Archive & ar, gurls::GurlsOptionsList& opt, const unsigned int /* file_version */)
{
    gurls::OptTypes type;
    std::string name;

    ar & type;
    ar & name;

    opt.setName(name);
    opt.removeOpt("Name");

    int n = 0;
    ar & n;

    for (int i = 0; i<n; ++i)
    {
        ar & name;
        ar & type;

        switch(type)
        {

        case gurls::StringOption:
        {
            gurls::OptString* tmp = new gurls::OptString();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::StringListOption:
        {
            gurls::OptStringList* tmp = new gurls::OptStringList();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::NumberOption:
        {
            gurls::OptNumber* tmp = new gurls::OptNumber();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::NumberListOption:
        {
            gurls::OptNumberList* tmp = new gurls::OptNumberList();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::FunctionOption:
        {
            gurls::OptFunction* tmp = new gurls::OptFunction();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::MatrixOption:
        {

            gurls::OptMatrixBase::MatrixType matType;
            ar & matType;

            gurls::OptMatrixBase* tmp;

#ifdef _BGURLS

            bool bigArray;
            ar & bigArray;

            if(bigArray)
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        tmp = new gurls::OptMatrix<gurls::BigArray<float> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<float> >*>(tmp));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        tmp = new gurls::OptMatrix<gurls::BigArray<double> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<double> >*>(tmp));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        tmp = new gurls::OptMatrix<gurls::BigArray<unsigned long> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::BigArray<unsigned long> >*>(tmp));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }
            else
#endif
            {
                switch(matType)
                {
                    case gurls::OptMatrixBase::FLOAT:
                        tmp = new gurls::OptMatrix<gurls::gMat2D<float> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<float> >*>(tmp));
                        break;
                    case gurls::OptMatrixBase::DOUBLE:
                        tmp = new gurls::OptMatrix<gurls::gMat2D<double> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<double> >*>(tmp));
                        break;
                    case gurls::OptMatrixBase::ULONG:
                        tmp = new gurls::OptMatrix<gurls::gMat2D<unsigned long> >();
                        ar & *(static_cast<gurls::OptMatrix<gurls::gMat2D<unsigned long> >*>(tmp));
                        break;
                    default:
                        throw gurls::gException(gurls::Exception_Unsupported_MatrixType);
                }
            }

            opt.addOpt(name, tmp);
            break;

        }
        case gurls::TaskSequenceOption:
        {
            gurls::OptTaskSequence* tmp = new gurls::OptTaskSequence();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::ProcessOption:
        {
            gurls::OptProcess* tmp = new gurls::OptProcess();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::OptListOption:
        {
            gurls::GurlsOptionsList* tmp = new gurls::GurlsOptionsList("", false);
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        case gurls::OptArrayOption:
        {
            gurls::OptArray* tmp = new gurls::OptArray();
            ar & (*tmp);
            opt.addOpt(name, tmp);
            break;
        }
        default:
            throw gurls::gException(gurls::Exception_Unknown_Option);

        }
    }

    opt.setName(opt.getName());
}



/**
  * Serializes an OptFunction to a generic archive
  */
template<class Archive>
void save(Archive & ar, const gurls::OptFunction& opt, const unsigned int /* file_version */)
{
    gurls::OptTypes type = opt.getType();
    ar & type;

    std::string name = opt.getName();
    ar & name;
}

/**
  * Deserializes an OptFunction from a generic archive
  */
template<class Archive>
void load(Archive & ar, gurls::OptFunction& opt, const unsigned int /* file_version */)
{
    gurls::OptTypes type;
    ar & type;

    std::string name;
    ar & name;

    opt = gurls::OptFunction(name);
}



template<class Archive, class OptType>
inline void serialize_opt( Archive& ar, OptType& opt, const unsigned int /*file_version*/)
{
    gurls::OptTypes type = opt.getType();
    ar & type;

    typename OptType::ValueType& value = opt.getValue();
    ar & value;
}

template<class Archive, class OptType>
void save_vector_opt(Archive & ar, const OptType& opt, const unsigned int /* file_version */)
{
    gurls::OptTypes type = opt.getType();
    ar & type;

    const typename OptType::ValueType& value = opt.getValue();

    int n = value.size();
    ar & n;

    for (int i = 0; i < n; ++i)
        ar & value[i];
}

template<class Archive, class OptType>
void load_vector_opt(Archive & ar, OptType& opt, const unsigned int /* file_version */)
{
    gurls::OptTypes type = opt.getType();
    ar & type;

    int n = 0;
    ar & n;

    opt.clear();

    typedef typename OptType::ValueType ValueType;
    typename ValueType::value_type tmp;

    for(int i = 0; i < n; ++i)
    {
        ar & tmp;
        opt << tmp;
    }

//    typedef typename OptType::ValueType ValueType;
//    ValueType value;

//    typename ValueType::value_type tmp;

//    for(int i = 0; i < n; ++i)
//    {
//        ar & tmp;
//        value.push_back(tmp);
//    }

//    opt.setValue(value);
}


/**
  * Serializes and deserializes an OptMatrix to/from a generic archive
  */
template<class Archive, typename T>
inline void serialize(Archive& ar, gurls::OptMatrix<gurls::gMat2D<T> >& opt, const unsigned int file_version)
{
    serialize_opt<Archive, gurls::OptMatrix<gurls::gMat2D<T> > >(ar, opt, file_version);
}

#ifdef _BGURLS
/**
  * Serializes and deserializes an OptMatrix to/from a generic archive
  */
template<class Archive, typename T>
inline void serialize(Archive& ar, gurls::OptMatrix<gurls::BigArray<T> >& opt, const unsigned int file_version)
{
    serialize_opt<Archive, gurls::OptMatrix<gurls::BigArray<T> > >(ar, opt, file_version);
}
#endif

/**
  * Serializes and deserializes an OptString to/from a generic archive
  */
template<class Archive>
void serialize(Archive& ar, gurls::OptString& opt, const unsigned int file_version)
{
    serialize_opt<Archive, gurls::OptString>(ar, opt, file_version);
}



/**
  * Serializes an OptStringList to a generic archive
  */
template<class Archive>
void save(Archive & ar, const gurls::OptStringList& opt, const unsigned int file_version)
{
    save_vector_opt<Archive, gurls::OptStringList>(ar, opt, file_version);
}

/**
  * Deserializes an OptStringList from a generic archive
  */
template<class Archive>
void load(Archive & ar, gurls::OptStringList& opt, const unsigned int file_version )
{
    load_vector_opt<Archive, gurls::OptStringList>(ar, opt, file_version);
}


/**
  * Serializes and deserializes an OptNumber to/from a generic archive
  */
template<class Archive>
void serialize(Archive & ar, gurls::OptNumber& opt, const unsigned int file_version)
{
    serialize_opt<Archive, gurls::OptNumber>(ar, opt, file_version);
}


/**
  * Serializes an OptNumberList to a generic archive
  */
template<class Archive>
void save(Archive & ar, const gurls::OptNumberList& opt, const unsigned int file_version)
{
    save_vector_opt<Archive, gurls::OptNumberList>(ar, opt, file_version);
}

/**
  * Deserializes an OptNumberList from a generic archive
  */
template<class Archive>
void load(Archive & ar, gurls::OptNumberList& opt, const unsigned int file_version)
{
    load_vector_opt<Archive, gurls::OptNumberList>(ar, opt, file_version);
}


/**
  * Serializes an OptTaskSequence to a generic archive
  */
template<class Archive>
void save(Archive & ar, const gurls::OptTaskSequence& opt, const unsigned int file_version)
{
    save_vector_opt<Archive, gurls::OptTaskSequence>(ar, opt, file_version);
}

/**
  * Deserializes an OptTaskSequence from a generic archive
  */
template<class Archive>
void load(Archive & ar, gurls::OptTaskSequence& opt, const unsigned int file_version)
{
    load_vector_opt<Archive, gurls::OptTaskSequence>(ar, opt, file_version);
}


/**
  * Serializes an OptProcess to a generic archive
  */
template<class Archive>
void save(Archive & ar, const gurls::OptProcess& opt, const unsigned int file_version)
{
    save_vector_opt<Archive, gurls::OptProcess>(ar, opt, file_version);
}

/**
  * Deserializes an OptProcess from a generic archive
  */
template<class Archive>
void load(Archive & ar, gurls::OptProcess& opt, const unsigned int file_version)
{
    load_vector_opt<Archive, gurls::OptProcess>(ar, opt, file_version);
}

} // namespace serialization
} // namespace boost

BOOST_SERIALIZATION_SPLIT_FREE(gurls::OptArray)
BOOST_SERIALIZATION_SPLIT_FREE(gurls::GurlsOptionsList)
BOOST_SERIALIZATION_SPLIT_FREE(gurls::OptFunction)
BOOST_SERIALIZATION_SPLIT_FREE(gurls::OptStringList)
BOOST_SERIALIZATION_SPLIT_FREE(gurls::OptNumberList)
BOOST_SERIALIZATION_SPLIT_FREE(gurls::OptTaskSequence)
BOOST_SERIALIZATION_SPLIT_FREE(gurls::OptProcess)


#endif // _GURLS_OPTSERIALIZATION_H_
