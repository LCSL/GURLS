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

#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>

#include "options.h"
#include "optlist.h"

using namespace std;

namespace gurls{

void GurlsOptionsList::setName(std::string newname){

    this->name = newname;
    this->table->erase("Name");
    this->table->insert(pair<std::string, GurlsOption*> ("Name", new OptString(newname)));
}

GurlsOptionsList::GurlsOptionsList(std::string ExpName, bool usedefopt): GurlsOption(OptListOption), name(ExpName){

    string key ("Name");
    GurlsOption* value;
    value = new OptString(ExpName);
	table = new std::map<std::string, GurlsOption* >();
    table->insert( pair<std::string,GurlsOption*>(key, value) );

    if(usedefopt){

        //		opt.combineclasses = @mean; % How to combine performance measure per class (mean/median/min/max?)
        //this->table["combineclasses"] = new OptFunction("mean", mean);
        (*table)["combineclasses"] = new OptFunction("mean");

        (*table)["name"] = new OptString(ExpName);
        (*table)["plotstr"] = new OptString(ExpName);

        //		WARNING: THE FILE EXTENSION AND THE POSSIBILITY TO SAVE
        //		THE EXPERIMENTS HAVE NOT BEEN IMPLEMENTED YET
#ifdef USE_BINARY_ARCHIVES
        (*table)["savefile"] = new OptString(ExpName.append(".bin"));
#else
        (*table)["savefile"] = new OptString(ExpName.append(".txt"));
#endif

        // ================================================== Algorithm options

        //		opt.kernel.type = 'rbf';
        //(*table)["singlelambda"] = new OptFunction("median", median);
        //(*table)["singlelambda"] = new OptFunction("median", median);
        (*table)["singlelambda"] = new OptFunction("median");
        (*table)["predbagmethod"] = new OptString("vote");
        // NOTE: lambda is searched between
        // [min(eig_r, opt.smallnumber), eig_1],
        // where r = rank, eig_1 = max eig val.
        (*table)["smallnumber"] = new OptNumber(1e-8);

        // ================================================== Directory options
        (*table)["tmpdir"] = new OptString(ExpName);

        // ===================================================== Output options
        (*table)["savekernel"] = new OptNumber(1);
        (*table)["saveanalysis"] = new OptNumber(1);
        //		opt.hoperf = @perf_precrec;
        (*table)["ploteval"] = new OptString("acc");
        //		WARNING: this should be an array of strings...
        (*table)["perfeval"] = new OptString("acc");

        // ======================================================== Data option
        (*table)["nholdouts"] = new OptNumber(1);
        (*table)["hoproportion"] = new OptNumber(0.2);
//        (*table)["nlambda"] = new OptNumber(100);
//        (*table)["nsigma"] =  new OptNumber(25);
        (*table)["nlambda"] = new OptNumber(20);
        (*table)["nsigma"] =  new OptNumber(10);


    // ======================================================== Pegasos option
        (*table)["subsize"]   = new OptNumber(50);
        (*table)["calibfile"] = new OptString("pippo");
        (*table)["epochs"]   = new OptNumber(4);

        // ============================================================== Quiet
        // Currenty either 0 or 1; levels of verbosity may be implemented later;
        (*table)["verbose"] = new OptNumber(1);

        // ======================================================= Version info
        (*table)["version"] = new OptString("0.1");

    }

}

GurlsOptionsList::~GurlsOptionsList()
{
    std::map<std::string, GurlsOption* >::iterator it, end;
    for(it = table->begin(), end = table->end(); it != end; ++it)
        delete (it->second);

    table->clear();
	delete table;
}

void GurlsOptionsList::printAll(){
    std::cout << *this;
}

void GurlsOptionsList::removeOpt(string key, bool deleteMembers)
{
    if(hasOpt(key))
    {
        if (deleteMembers)
            delete (*table)[key];

        table->erase(key);
    }
}
std::ostream& GurlsOptionsList::operator<<(std::ostream& os){
    return os << *this;
}


GURLS_EXPORT std::ostream& operator<<(std::ostream& os, GurlsOptionsList& opt) {
    std::map<std::string, GurlsOption* >::iterator it;
    GurlsOption* current;
    os << std::endl
       << "~~~~~~~ GurlsOptionList: " << opt.getName()
       << std::endl;
    for (it = opt.table->begin(); it != opt.table->end(); it++){
        current = (*it).second;
        os << "\t[ " << (*it).first << " ] = ";
        current->operator <<(os);
        os << endl;
    }
    os << "~~~~~~~";
    return os;
}

bool GurlsOptionsList::addOpt(std::string key, GurlsOption* value){
    table->insert( pair<std::string,GurlsOption*>(key, value) );
    return true;
}

bool GurlsOptionsList::addOpt(std::string key, std::string value){
    OptString* v = new OptString(value);
    table->insert( pair<std::string,GurlsOption*>(key, v) );
    return true;
}

GurlsOption* GurlsOptionsList::getOpt(std::string key){
    GurlsOption* gout = (*table)[key];
    if (!gout){
        throw gException(Exception_Parameter_Not_Definied_Yet + "( " + key + " )");
    }
    return gout;
}

std::string GurlsOptionsList::getOptAsString(std::string key){

    GurlsOption* opt = this->getOpt(key);
    if (opt->getType() == StringOption){
        return static_cast<OptString*>(opt)->getValue();
    }else {
        throw gException(Exception_Unknown_Option);
    }
}


double GurlsOptionsList::getOptAsNumber(std::string key){

    GurlsOption* opt = this->getOpt(key);
    if (opt->getType() == NumberOption){
        return static_cast<OptNumber*>(opt)->getValue();
    }else {
        throw gException(Exception_Unknown_Option);
    }
}

void GurlsOptionsList::save(const std::string& fileName) const
{
    std::ofstream outstream(fileName.c_str());

    if(!outstream.is_open())
        throw gException("Could not open file " + fileName);

    oarchive outar(outstream);
    outar << *this;

    outstream.close();
}

void GurlsOptionsList::load(const std::string& fileName)
{
    std::ifstream instream(fileName.c_str());

    if(!instream.is_open())
        throw gException("Could not open file " + fileName);

    iarchive inar(instream);
    inar >> *this;

    instream.close();
}


}

