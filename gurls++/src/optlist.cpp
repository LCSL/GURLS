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
    this->table.erase("Name");
    this->table.insert(pair<std::string, GurlsOption*> ("Name", new OptString(newname)));
}

GurlsOptionsList::GurlsOptionsList(std::string ExpName, bool usedefopt): GurlsOption(OptListOption), name(ExpName){

    string key ("Name");
    GurlsOption* value;
    value = new OptString(ExpName);
    this->table.insert( pair<std::string,GurlsOption*>(key, value) );

    if(usedefopt){

        //		opt.combineclasses = @mean; % How to combine performance measure per class (mean/median/min/max?)
        //this->table["combineclasses"] = new OptFunction("mean", mean);
        this->table["combineclasses"] = new OptFunction("mean");

        this->table["plotstr"] = new OptString(ExpName);

        //		WARNING: THE FILE EXTENSION AND THE POSSIBILITY TO SAVE
        //		THE EXPERIMENTS HAVE NOT BEEN IMPLEMENTED YET
        this->table["savefile"] = new OptString(ExpName.append(".mat"));

        // ================================================== Algorithm options

        //		opt.kernel.type = 'rbf';
        //this->table["singlelambda"] = new OptFunction("median", median);
        //this->table["singlelambda"] = new OptFunction("median", median);
        this->table["singlelambda"] = new OptFunction("median");
        this->table["predbagmethod"] = new OptString("vote");
        // NOTE: lambda is searched between
        // [min(eig_r, opt.smallnumber), eig_1],
        // where r = rank, eig_1 = max eig val.
        this->table["smallnumber"] = new OptNumber(1e-8);

        // ================================================== Directory options
        this->table["tmpdir"] = new OptString(ExpName);

        // ===================================================== Output options
        this->table["savekernel"] = new OptNumber(1);
        this->table["saveanalysis"] = new OptNumber(1);
        //		opt.hoperf = @perf_precrec;
        this->table["ploteval"] = new OptString("acc");
        //		WARNING: this should be an array of strings...
        this->table["perfeval"] = new OptString("acc");

        // ======================================================== Data option
        this->table["nholdouts"] = new OptNumber(1);
        this->table["hoproportion"] = new OptNumber(0.2);
        this->table["nlambda"] = new OptNumber(100);
        this->table["nsigma"] =  new OptNumber(25);

        // ============================================================== Quiet
        // Currenty either 0 or 1; levels of verbosity may be implemented later;
        this->table["verbose"] = new OptNumber(1);

        // ======================================================= Version info
        this->table["version"] = new OptString("0.1");

    }

}

void GurlsOptionsList::printAll(){
    std::cout << *this;
}

void GurlsOptionsList::removeOpt(string key, bool deleteMembers)
{
    if(hasOpt(key))
    {
        if (deleteMembers)
            delete table[key];

        table.erase(key);
    }
}
std::ostream& GurlsOptionsList::operator<<(std::ostream& os){
    return os << *this;
}


std::ostream& operator<<(std::ostream& os, GurlsOptionsList& opt) {
    std::map<std::string, GurlsOption* >::iterator it;
    GurlsOption* current;
    os << std::endl
       << "~~~~~~~ GurlsOptionList: " << opt.getName()
       << std::endl;
    for (it = opt.table.begin(); it != opt.table.end(); it++){
        current = (*it).second;
        os << "\t[ " << (*it).first << " ] = ";
        current->operator <<(os);
        os << endl;
    }
    os << "~~~~~~~";
    return os;
}

bool GurlsOptionsList::addOpt(std::string key, GurlsOption* value){
    this->table.insert( pair<std::string,GurlsOption*>(key, value) );
}

bool GurlsOptionsList::addOpt(std::string key, std::string value){
    OptString* v = new OptString(value);
    this->table.insert( pair<std::string,GurlsOption*>(key, v) );
}

GurlsOption* GurlsOptionsList::getOpt(std::string key){
    GurlsOption* gout = this->table[key];
    if (!gout){
        throw gException(Exception_Parameter_Not_Definied_Yet);
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


}

