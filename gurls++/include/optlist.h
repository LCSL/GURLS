/*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011, IIT@MIT Lab
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

#ifndef _GURLS_OPTLIST_H_
#define _GURLS_OPTLIST_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <typeinfo>

#include "exports.h"
#include "gmat2d.h"
#include "gvec.h"
#include "options.h"
#include "exceptions.h"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>


namespace gurls {



#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4251)
#endif

/**
  * \ingroup Settings
  * \brief GurlsOptionsList is an option containing a list of options
  * mapped by name.
  */
class GURLS_EXPORT GurlsOptionsList: public GurlsOption
{
private:
    std::string name;                               ///< Option name
    std::map<std::string, GurlsOption* >* table;    ///< Options list, indexed by name

public:


    /**
      * Constructor. Builds an optionlist with a name and optionally a set of default options
      *
      * \param ExpName name of the options list
      * \param usedefopt if \a true the list is filled with a set of default options, if \a false the list is left empty
      */
    GurlsOptionsList(std::string ExpName, bool usedefopt = false);


    /**
      * Destructor
      */
    ~GurlsOptionsList();


    /**
      * Adds a generic option to the list indexed with a specified key
      */
    bool addOpt(std::string key, GurlsOption* value);

    /**
      * Adds a string option to the list indexed with a specified key
      */
    bool addOpt(std::string key, std::string value);

    /**
      * Returns a pointer to a generic option mapped with a key
      */
    GurlsOption* getOpt(std::string key);

    /**
      * Returns a pointer to a generic option mapped with a key
      */
    const GurlsOption* getOpt(std::string key) const;

    /**
      * Returns a pointer to a T option mapped with a key
      */
    template<class T>
    T* getOptAs(std::string key)
    {
        return T::dynacast(this->getOpt(key));
    }

    /**
      * Returns a pointer to a T option mapped with a key
      */
    template<class T>
    const T* getOptAs(std::string key) const
    {
        return T::dynacast(this->getOpt(key));
    }

    /**
      * Returns a reference to the value contained into an option mapped with a key
      */
    template<class T>
    typename T::ValueType& getOptValue(std::string key)
    {
        return this->getOptAs<T>(key)->getValue();
    }

    /**
      * Returns a reference to the value contained into an option mapped with a key
      */
    template<class T>
    const typename T::ValueType& getOptValue(std::string key) const
    {
        return this->getOptAs<T>(key)->getValue();
    }

    /**
      * Returns a string option mapped with a key
      */
    std::string getOptAsString(std::string key) const;

    /**
      * Returns the list name
      */
    std::string getName() const {return this->name;}

    /**
      * Sets the list name
      */
    void setName(std::string);

    /**
      * Returns a numeric option mapped with a key
      */
    double getOptAsNumber(std::string key) const;

    /**
      * Prints the options list
      */
    void printAll();

    /**
      * Checks if the list has an option mapped with a specified key
      */
    bool hasOpt(std::string key) const {return table->count(key)>0;}

    /**
      * Removes the option mapped with a specified key
      *
      * \param key string key
      * \param deleteMembers If true deallocates the removed option, if false option will be only detached
      */
    void removeOpt(std::string key, bool deleteMembers = true);

    /**
      * Checks if the option has the given type
      */
    virtual bool isA(OptTypes id) const { return (id == OptListOption); }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref GurlsOptionsList
      */
    static GurlsOptionsList* dynacast(GurlsOption* opt) {
        if (opt->isA(OptListOption) ){
            return static_cast<GurlsOptionsList*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Tries to cast a pointer to a generic option to a pointer to an \ref GurlsOptionsList
      */
    static const GurlsOptionsList* dynacast(const GurlsOption* opt)
    {
        if (opt->isA(OptListOption) ){
            return static_cast<const GurlsOptionsList*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }

    /**
      * Returns the number of options on the list
      */
    int size() const {return table->size();}

    /**
      * Returns a pointer to the idx-th option into the list
      */
    GurlsOption* operator[] (int idx){

        if ( idx > this->size() ) {
            throw gException(gurls::Exception_Index_Out_of_Bound);
        }

        std::map<std::string, GurlsOption* >::iterator itr = table->begin();
        std::map<std::string, GurlsOption* >::iterator end = table->end();
        for (int i = 0; itr!=end; ++i, ++itr){
            /* Do nothing else then following the iterator */
        }
        return itr->second;
    }

    /**
      * Writes a GurlsOptionsList to a stream
      */
    friend GURLS_EXPORT std::ostream& operator<<(std::ostream& os, GurlsOptionsList& opt);

    /**
      * Writes the list to a stream
      */
    virtual std::ostream& operator<<(std::ostream& os);

    /**
      *
      */
    template<typename T>
    void copyOpt(std::string key, const GurlsOptionsList &from)
    {
        const GurlsOption* toCopy = from.getOpt(key);

        GurlsOption* newOpt = NULL;

        switch(toCopy->getType())
        {
        case StringOption:
            newOpt = new OptString(OptString::dynacast(toCopy)->getValue());
            break;
        case NumberOption:
            newOpt = new OptNumber(OptNumber::dynacast(toCopy)->getValue());
            break;
        case StringListOption:
            newOpt = new OptStringList(OptStringList::dynacast(toCopy)->getValue());
            break;
        case NumberListOption:
            newOpt = new OptNumberList(OptNumberList::dynacast(toCopy)->getValue());
            break;
        case FunctionOption:
            newOpt = new OptFunction(OptFunction::dynacast(toCopy)->getName());
            break;
        case MatrixOption:
        case VectorOption:
        {
           const OptMatrixBase* base = dynamic_cast<const OptMatrixBase*>(toCopy);

            if(base == NULL)
                throw gException(Exception_Illegal_Dynamic_Cast);

            if(base->getMatrixType() == OptMatrixBase::ULONG)
            {
                const gMat2D<unsigned long> & mat = OptMatrix<gMat2D<unsigned long> >::dynacast(toCopy)->getValue();

                gMat2D<unsigned long>* newMat = new gMat2D<unsigned long>(mat);
                newOpt = new OptMatrix<gMat2D<unsigned long> >(*newMat);
            }
            else
            {
                const gMat2D<T> & mat = OptMatrix<gMat2D<T> >::dynacast(toCopy)->getValue();
                gMat2D<T>* newMat = new gMat2D<T>(mat);
                newOpt = new OptMatrix<gMat2D<T> >(*newMat);
            }

        }
            break;
        case OptListOption:
        {
            const GurlsOptionsList* toCopy_list = GurlsOptionsList::dynacast(toCopy);

            GurlsOptionsList* list = new GurlsOptionsList(toCopy_list->name);

            std::map<std::string, GurlsOption* >::const_iterator it, end;

            list->removeOpt("Name");

            for(it = toCopy_list->table->begin(), end = toCopy_list->table->end(); it != end; ++it)
                list->copyOpt<T>(it->first, *toCopy_list);

            list->setName(list->getOptAsString("Name"));

            newOpt = list;
        }
            break;
        case TaskSequenceOption:
            newOpt = new OptTaskSequence(OptTaskSequence::dynacast(toCopy)->getValue());
            break;
        case TaskIDOption:
        case GenericOption:
            break;
        }

        if(newOpt != NULL)
            addOpt(key, newOpt);
    }


    /**
      * Serializes the list to file
      */
    void save(const std::string& fileName) const;

    /**
      * Deserializes the list from file
      */
    void load(const std::string& fileName);

    friend class boost::serialization::access;

    /**
      * Serializes the option to a generic archive
      */
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        ar & this->name;
        int n = table->size();
        int type = -1;
        ar & n;
        std::map<std::string, GurlsOption* >::const_iterator itr = table->begin();
        std::map<std::string, GurlsOption* >::const_iterator end = table->end();
        for (int i = 0; itr!=end; ++i, ++itr){
            std::string s = itr->first;
            ar & s;
            GurlsOption* opt0 = itr->second;
            type = opt0->getType();
            ar & type;
//			std::cout << " ° " << s << ": " << type << std::endl;
            if (type == StringOption){
                OptString* opt = static_cast<OptString*>(opt0);
                ar & (*opt);
            } else if (type == StringListOption){
                OptStringList* opt = static_cast<OptStringList*>(opt0);
                ar & (*opt);
            } else if (type == NumberOption){
                OptNumber* opt = static_cast<OptNumber*>(opt0);
                ar & (*opt);
            } else if (type == NumberListOption){
                OptNumberList* opt = static_cast<OptNumberList*>(opt0);
                ar & (*opt);
            } else if (type == FunctionOption){
                OptFunction* opt = static_cast<OptFunction*>(opt0);
                ar & (*opt);
            } else if (type == MatrixOption) {

                OptMatrixBase* optM = static_cast<OptMatrixBase*>(opt0);
                OptMatrixBase::MatrixType matType = optM->getMatrixType();

                ar & matType;

                switch(matType)
                {
                    case OptMatrixBase::FLOAT:
                        ar & *(static_cast<OptMatrix<gMat2D<float> >*>(optM));
                        break;
                    case OptMatrixBase::DOUBLE:
                        ar & *(static_cast<OptMatrix<gMat2D<double> >*>(optM));
                        break;
                    case OptMatrixBase::ULONG:
                        ar & *(static_cast<OptMatrix<gMat2D<unsigned long> >*>(optM));
                        break;
                    default:
                        throw gException(Exception_Unsupported_MatrixType);
                }

            } else if (type == TaskSequenceOption) {
                OptTaskSequence* opt = static_cast<OptTaskSequence*>(opt0);
                ar & (*opt);
            } else if (type == OptListOption){
                GurlsOptionsList* opt = GurlsOptionsList::dynacast(opt0);
                ar & (*opt);
            } else {
                // AN EXCEPTION SHOULD BE RAISED
            }


            //ar & (*opt);
//			std::cout << " done." << std::endl;
        }
    }

    /**
      * Deserializes the option from a generic archive
      */
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        ar & this->name;
        removeOpt("Name");
        int n = 0;
        int type = -1;
        ar & n;
        //GurlsOption* opt;
        for (int i = 0; i<n; ++i){
            //GurlsOption* opt;
            std::string s;
            ar & s;
            ar & type;
            if (type == StringOption){
                OptString* opt = new OptString();
                ar & (*opt);
                this->addOpt(s, opt);
            } else if (type == StringListOption){
                OptStringList* opt = new OptStringList();
                ar & (*opt);
                this->addOpt(s, opt);
            } else if (type == NumberOption){
                OptNumber* opt = new OptNumber();
                ar & (*opt);
                this->addOpt(s, opt);
            } else if (type == NumberListOption){
                OptNumberList* opt = new OptNumberList();
                ar & (*opt);
                this->addOpt(s, opt);
            } else if (type == FunctionOption){
                OptFunction* opt = new OptFunction("");
                ar & (*opt);
                this->addOpt(s, opt);
            } else if (type == MatrixOption) {

                OptMatrixBase::MatrixType matType;
                ar & matType;

                OptMatrixBase* opt;

                switch(matType)
                {
                    case OptMatrixBase::FLOAT:
                        opt = new OptMatrix<gMat2D<float> >();
                        ar & *(static_cast<OptMatrix<gMat2D<float> >*>(opt));
                        break;
                    case OptMatrixBase::DOUBLE:
                        opt = new OptMatrix<gMat2D<double> >();
                        ar & *(static_cast<OptMatrix<gMat2D<double> >*>(opt));
                        break;
                    case OptMatrixBase::ULONG:
                        opt = new OptMatrix<gMat2D<unsigned long> >();
                        ar & *(static_cast<OptMatrix<gMat2D<unsigned long> >*>(opt));
                        break;
                    default:
                        throw gException(Exception_Unsupported_MatrixType);
                }

                this->addOpt(s, opt);

            } else if (type == TaskSequenceOption) {
                OptTaskSequence* opt = new OptTaskSequence();
                ar & (*opt);
                this->addOpt(s, opt);
            } else if (type == OptListOption){
                GurlsOptionsList* opt = new GurlsOptionsList("tmp", false);
                ar & (*opt);
                this->addOpt(s, opt);
            } else {
                // AN EXCEPTION SHOULD BE RAISED
            }

            //ar & (*opt);
            //this->addOpt(s, opt);
        }
        this->setName(this->name);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()


};

#ifdef _WIN32
#pragma warning(pop)
#endif

}

#endif // _GURLS_OPTLIST_H_
