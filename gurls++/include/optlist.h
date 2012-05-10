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

#include "gmat2d.h"
#include "gvec.h"
#include "options.h"
#include "exceptions.h"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>


namespace gurls {

class GurlsOptionsList: public GurlsOption
{
private:
    std::string name;
    std::map<std::string, GurlsOption* > table;

public:
    GurlsOptionsList(std::string ExpName, bool usedefopt = false);
    bool addOpt(std::string key, GurlsOption* value);
    bool addOpt(std::string key, std::string value);
    GurlsOption* getOpt(std::string key);
    std::string getOptAsString(std::string key);
    std::string getName() const {return this->name;}
    void setName(std::string);
    double getOptAsNumber(std::string key);
    void printAll();
    bool hasOpt(std::string key) {return table.count(key)>0;}
    void removeOpt(std::string key, bool deleteMembers = true);


    virtual bool isA(int id) { return (id == OptListOption); }

    static GurlsOptionsList* dynacast(GurlsOption* opt) {
        if (opt->isA(OptListOption) ){
            return static_cast<GurlsOptionsList*>(opt);
        } else {
            throw gException(gurls::Exception_Illegal_Dynamic_Cast);
        }
    }
    int size() const {return table.size();}

    GurlsOption* operator[] (int idx){

        if ( idx > this->size() ) {
            throw gException(gurls::Exception_Index_Out_of_Bound);
        }

        std::map<std::string, GurlsOption* >::iterator itr = table.begin();
        std::map<std::string, GurlsOption* >::iterator end = table.end();
        for (int i = 0; i<=idx, itr!=end; i++, itr++){
            /* Do nothing else then following the iterator */
        }
        return itr->second;
    }

    friend std::ostream& operator<<(std::ostream& os, GurlsOptionsList& opt);
    virtual std::ostream& operator<<(std::ostream& os);


    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const{
        ar & this->type;
        ar & this->name;
        int n = table.size();
        int type = -1;
        ar & n;
        std::map<std::string, GurlsOption* >::const_iterator itr = table.begin();
        std::map<std::string, GurlsOption* >::const_iterator end = table.end();
        for (int i = 0; i<n, itr!=end; i++, itr++){
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
                OptMatrix<gMat2D<float> >* opt = static_cast<OptMatrix<gMat2D<float> >*>(opt0);
                ar & (*opt);
            } else if (type == TaskSequenceOption) {
                OptTaskSequence* opt = static_cast<OptTaskSequence*>(opt0);
                ar & (*opt);
            } else if (type == OptListOption){
                GurlsOptionsList* opt = static_cast<GurlsOptionsList*>(opt0);
                ar & (*opt);
            } else {
                // AN EXCEPTION SHOULD BE RAISED
            }


            //ar & (*opt);
//			std::cout << " done." << std::endl;
        }
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar & this->type;
        ar & this->name;
        this->setName(this->name);
        int n = 0;
        int type = -1;
        ar & n;
        //GurlsOption* opt;
        for (int i = 0; i<n; i++){
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
                OptMatrix<gMat2D<float> >* opt = new OptMatrix<gMat2D<float> >();
                ar & (*opt);
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
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()


};

}

#endif // _GURLS_OPTLIST_H_
