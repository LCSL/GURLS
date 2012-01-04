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


#ifndef _GURLS_EXCEPTIONS_H_
#define _GURLS_EXCEPTIONS_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <exception>

namespace gurls {

/*
  STATIC GLOBAL FIELDS DEFINED IN THE NAMESPACE TO IDENTIFY TYPICAL EXCEPTIONS
  */

static std::string
Exception_Incipit("ERROR! ");

static std::string
Exception_Wrong_Ownership(Exception_Incipit+"An attempt to modify an array that is not the owner of its data occurred. ");

static std::string
Exception_Wrong_Memory_Access(Exception_Incipit+"An attempt to acces a non-existent memory location occurred. ");

static std::string
Exception_Inconsistent_Size(Exception_Incipit+"An attempt to combine arrays with inconsistent dimensions occurred. ");

static std::string
Exception_Functionality_Not_Implemented(Exception_Incipit+"An attempt to use a functionality that is not implemented yet occurred. ");

static std::string
Exception_Square_Matrix_Required(Exception_Incipit+"An attempt to use a general matrix instead of the required square matrix occurred. ");

static std::string
Exception_Unknown_Option(Exception_Incipit+"An unknown option has been used.");

static std::string
Exception_Required_Parameter_Missing(Exception_Incipit+"One of the parameters required to run the algorithm is missing.");

static std::string
Exception_Parameter_Not_Definied_Yet(Exception_Incipit+"The requested parameter has not been defined yet.");

static std::string
Exception_Logical_Operator(Exception_Incipit+"An unknown logical comparison has been required.");

static std::string
Exception_Illegal_Dynamic_Cast(Exception_Incipit+"An illegal dynamic cast occured.");

static std::string
Exception_Illegat_Argument_Value(Exception_Incipit+"The value of the input variable is not allowed.");

static std::string
Exception_Invalid_Reshape_Arguments(Exception_Incipit+"To RESHAPE the number of elements must not change.");

static std::string
Exception_Index_Out_of_Bound(Exception_Incipit+"Index exceeds matrix dimensions.");

static std::string
Exception_Gurls_Inconsistent_Processes_Number(Exception_Incipit+"The number of elements in the list of processes/tasks is not consistent. ");

static std::string
Exception_Gurls_Invalid_ProcessID(Exception_Incipit+"Invalid process ID. ");

static std::string
Exception_Invalid_TaskSequence(Exception_Incipit+"Invalid task name specification. ");

static std::string
Exception_Unknown_Function(Exception_Incipit+"Unknown function.");



/*
  Class designed to deal with exceptions in Gurls++ package.
  The current mechanism is very simplistic: all the exceptions are instances of
  this base class and different error types are identified only by the message
  conveyed by that specific instance. Convenience static strings have been
  defined within the name space in order to identify typical and frequent
  exceptions (see above in this file).
  */
class gException: public std::exception {
private:
	std::string msg;
public:
	gException( std::string message ) : msg(message) { };
	virtual ~gException() throw () { };
	inline std::string getMessage() { return std::string("[")+this->what()+"]: " + msg; };
};

}

#endif // _GURLS_EXCEPTIONS_H_
