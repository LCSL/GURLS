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



#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <ctime>
#include <cassert>

#include "options.h"
#include "optlist.h"

#include <map>
#include<vector>
#include<list>


using namespace gurls;
using namespace std;

int main(int argc, char *argv[])
{

	GurlsOptionsList table("Test_Experiment", true);
	GurlsOptionsList nestedtable("NestedList", false);
	nestedtable.addOpt("pi greco", new OptNumber(3.14));

	string key("sigma");
	string value("1.56f");

	try{
		table.addOpt(key, value);
		string o2 = table.getOptAsString("Name");
		cout << o2 << endl;
		o2 = table.getOptAsString(key);
		cout << o2 << endl;

		string value2("verifica_funzionamento_inizializzazione_diretta_da_stringa");
		OptString opt2 = value2;
		table.addOpt("key_verifica", &opt2);
		cout << table.getOptAsString("key_verifica") << endl;
		opt2.setValue("magari_funziona_anche_questo?");
		cout << table.getOptAsString("key_verifica") << endl;



		double value3 = 2.3467;
		OptNumber opt3 = value3;
		table.addOpt("lambda", &opt3);
		cout << static_cast<OptNumber*>(table.getOpt("lambda"))->getValue() << endl;
		cout << table.getOptAsNumber("lambda") << endl;

		table.addOpt("nested", &nestedtable);

		//table.printAll();

		cout << table << endl;

	}
	catch (gException& e){
		cout << table.getOpt(key)->getType() << endl;
	}

	return 0;
}
