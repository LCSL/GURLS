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

#ifndef _GURLS_OPT_H_
#define _GURLS_OPT_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <typeinfo>

#include "gmat2d.h"
#include "gvec.h"
#include "exceptions.h"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>


namespace gurls {

enum OptTypes	{GenericOption, StringOption, NumberOption,
				 StringListOption, NumberListOption, FunctionOption,
				 MatrixOption, VectorOption,
				OptListOption, TaskSequenceOption, TaskIDOption};


/*
 * Typical functions used within the gurls package to combine or choose
 * among different regularization parameters or multiple numeric
 * options.
 */
double mean(double* v, int n);
double min(double* v, int n);
double max(double* v, int n);
double median(double* v, int n);


/*
	GURLSOPTION is an abstraction of a generic `option', which is
	widely used within the GURLS++ package to store either numeric
	parameters necessary to configure specific algorigms or sequences
	of strings holding the names of the specific procedures that
	have to be performed.

	The instances of the GURLSOPTION class hold information about
	the type (one of the elements in the OptTypes enumeration),
	while the value related to each specific option is stored using
	the attributes of the subclasses of GURLSOPTION.

  */
class GurlsOption
{
protected:
	int type;

public:
	GurlsOption(int t): type(t) {};
	const int getType() const {return type;}
	virtual ~GurlsOption(){};
	virtual bool isA(int id) { return (id == GenericOption); }
	const virtual std::type_info& getDataID(){
		return typeid(GurlsOption);
	}

	friend std::ostream& operator<<(std::ostream& os, GurlsOption& opt);
	virtual std::ostream& operator<<(std::ostream& os){return os;}

//	friend class boost::serialization::access;
//	template<class Archive>
//	void save(Archive & , const unsigned int) const;
//	template<class Archive>
//	void load(Archive & , const unsigned int) ;
//	BOOST_SERIALIZATION_SPLIT_MEMBER()

};

class OptString: public GurlsOption
{
private:
	std::string value;

public:
	OptString(): GurlsOption(StringOption), value(""){};
	OptString(const char* str): GurlsOption(StringOption),value(str){}
	OptString(std::string& str): GurlsOption(StringOption),value(str){}
	OptString& operator=(const OptString& other);

	virtual ~OptString(){ value.clear();}


	OptString& operator=(const std::string& other){
		this->type = StringOption;
		this->value = other;
		return *this;
	}

	void setValue(const std::string newvalue) {value = newvalue;}
	std::string& getValue() { return value;}

	virtual bool isA(int id) { return (id == StringOption); }
	static OptString* dynacast(GurlsOption* opt) {
		if (opt->isA(StringOption) ){
			return static_cast<OptString*>(opt);
		} else {
			throw gException(gurls::Exception_Illegal_Dynamic_Cast);
		}
	}

	virtual std::ostream& operator<<(std::ostream& os);

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int /* file_version */) const{
		ar & this->type;
		ar & this->value;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int /* file_version */){
		ar & this->type;
		ar & this->value;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};


class OptStringList: public GurlsOption
{
private:
	std::vector<std::string> value;

public:
	OptStringList(): GurlsOption(StringListOption), value(){};
	OptStringList(std::vector<std::string>& vec):
		GurlsOption(StringListOption), value(vec) {}
	OptStringList(std::string& str): GurlsOption(StringListOption),
		value(){ value.push_back(str); }
	OptStringList& operator=(const OptStringList& other);

	virtual ~OptStringList(){
		value.clear();
	}

	void setValue(const std::vector<std::string> newvalue) {
		value = newvalue;
	}
	void add(const std::string str){
		value.push_back(str);
	}

	std::vector<std::string>& getValue() { return value;}

	virtual bool isA(int id) { return (id == StringListOption); }
	static OptStringList* dynacast(GurlsOption* opt) {
		if (opt->isA(StringListOption) ){
			return static_cast<OptStringList*>(opt);
		} else {
			throw gException(gurls::Exception_Illegal_Dynamic_Cast);
		}
	}
	//friend std::ostream& operator<<(std::ostream& os, OptStringList& opt);
	virtual std::ostream& operator<<(std::ostream& os);

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int /* file_version */) const{
		ar & this->type;
		int n = this->value.size();
		ar & n;
		for (int i = 0; i < n; i++){
			ar & this->value.at(i);
		}
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int /* file_version */){
		ar & this->type;
		int n = 0;
		ar & n;
		value.clear();
		for(int i = 0; i < n; i++){
			std::string s("");
			ar & s;
			this->value.push_back(s);
		}
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()

};



class OptNumber: public GurlsOption
{
private:
	double value;

public:
	OptNumber(): GurlsOption(NumberOption), value(0) {}
	OptNumber(double v): GurlsOption(NumberOption), value(v) {}
	OptNumber& operator=(const OptNumber& other);

	virtual ~OptNumber(){};

	OptNumber& operator=(double other){
		this->type = NumberOption;
		this->value = other;
		return *this;
	}

	void setValue(double newvalue) {value = newvalue;}
	double& getValue() {return value;}

	virtual bool isA(int id) { return (id == NumberOption); }
	static OptNumber* dynacast(GurlsOption* opt) {
		if (opt->isA(NumberOption) ){
			return static_cast<OptNumber*>(opt);
		} else {
			throw gException(gurls::Exception_Illegal_Dynamic_Cast);
		}
	}
	// friend std::ostream& operator<<(std::ostream& os, OptNumber& opt);
	virtual std::ostream& operator<<(std::ostream& os);

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int /* file_version */) const{
		ar & this->type;
		ar & this->value;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int /* file_version */){
		ar & this->type;
		ar & this->value;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()

};

class OptNumberList: public GurlsOption
{
private:
	std::vector<double> value;

public:
	OptNumberList(): GurlsOption(NumberListOption), value(){};
	OptNumberList(std::vector<double>& vec):
		GurlsOption(NumberListOption), value(vec) {}
	OptNumberList(double v):
		GurlsOption(NumberListOption), value(){ value.push_back(v);}
	OptNumberList(double *v, int n):
		GurlsOption(NumberListOption), value(){
		int i = 0;
		while (i++<n){
			value.push_back(*v++);
		}
	}

	OptNumberList& operator=(const OptNumberList& other);

	virtual ~OptNumberList(){}

	void setValue(const std::vector<double> newvalue) {
		value = newvalue;
	}
	void add(const double d){
		value.push_back(d);
	}

	std::vector<double>& getValue() { return value;}

	virtual bool isA(int id) { return (id == NumberListOption); }
	static OptNumberList* dynacast(GurlsOption* opt) {
		if (opt->isA(NumberListOption) ){
			return static_cast<OptNumberList*>(opt);
		} else {
			throw gException(gurls::Exception_Illegal_Dynamic_Cast);
		}
	}
	// friend std::ostream& operator<<(std::ostream& os, OptNumberList& opt);
	virtual std::ostream& operator<<(std::ostream& os);

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int /* file_version */) const{
		ar & this->type;
		int n = this->value.size();
		ar & n;
		for (int i = 0; i < n; i++){
			ar & this->value.at(i);
		}
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int /* file_version */){
		ar & this->type;
		int n = 0;
		ar & n;
		value.clear();
		for(int i = 0; i < n; i++){
			double s;
			ar & s;
			this->value.push_back(s);
		}
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()

};


class OptFunction: public GurlsOption
{
private:
	//double (*value)(double*, int);
	std::string name;

public:
	OptFunction(std::string func_name): GurlsOption(FunctionOption), name(func_name) {}
	//OptFunction(std::string func_name): GurlsOption(FunctionOption), name(func_name), value(0) {}
	//OptFunction(std::string func_name, double (*v)(double*, int)): GurlsOption(FunctionOption), name(func_name), value(v) {}
	OptFunction& operator=(const OptFunction& other);

	virtual ~OptFunction(){};

//	void setValue(std::string func_name, double (*newvalue)(double*, int)) {
//		name = func_name;
//		value = newvalue;
//	}
	void setValue(std::string func_name) {
		name = func_name;
	}
	std::string getName() const {return name;}

	double& getValue(double* array, int n) {
		double *v = new double[1];
		//*v = (*value)(array,n);
		if (!name.compare("mean")){
			*v = mean(array, n);
		} else if(!name.compare("min")){
			*v = min(array, n);
		} else if(!name.compare("max")){
			*v = max(array, n);
		} else if(!name.compare("median")){
			*v = median(array, n);
		} else {
			*v = std::numeric_limits<double>::signaling_NaN();
			throw gException(gurls::Exception_Unknown_Function);
		}
		return *(v);
	}

	virtual bool isA(int id) { return (id == FunctionOption); }
	static OptFunction* dynacast(GurlsOption* opt) {
		if (opt->isA(FunctionOption) ){
			return static_cast<OptFunction*>(opt);
		} else {
			throw gException(gurls::Exception_Illegal_Dynamic_Cast);
		}
	}
	// friend std::ostream& operator<<(std::ostream& os, OptFunction& opt);
	virtual std::ostream& operator<<(std::ostream& os);

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int /* file_version */) const{
		ar & this->type;
		ar & this->name;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int /* file_version */){
		ar & this->type;
		ar & this->name;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};


template <typename Matrix>
//class OptMatrixTMP: public GurlsOption
class OptMatrix: public GurlsOption
{
private:
	Matrix& value;

public:
	OptMatrix(): GurlsOption(MatrixOption) , value (*(new Matrix(2,2))){}
	OptMatrix(Matrix& m): GurlsOption(MatrixOption), value(m) {}
	OptMatrix<Matrix>& operator=(const OptMatrix<Matrix>& other);

	virtual ~OptMatrix(){}

	OptMatrix& operator=(const Matrix& other){
		this->type = MatrixOption;
		this->value = other;
		return *this;
	}

	void setValue(const Matrix& newvalue) {value = newvalue;}
	Matrix& getValue() { return value;}

	virtual bool isA(int id) { return (id == MatrixOption); }
	static OptMatrix* dynacast(GurlsOption* opt) {
		if (opt->isA(MatrixOption) ){
			return static_cast<OptMatrix*>(opt);
		} else {
			throw gException(gurls::Exception_Illegal_Dynamic_Cast);
		}
	}

	//template <typename U>
	//friend std::ostream& operator<<(std::ostream& os, OptMatrixTMP<U>& opt);
	virtual std::ostream& operator<<(std::ostream& os);


	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int /* file_version */) const{
		ar & this->type;
		ar & this->value;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int /* file_version */){
		ar & this->type;
		ar & this->value;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()

};

static const std::string TASKDESC_SEPARATOR(":");

class OptTaskSequence: public GurlsOption
{
private:
	std::vector<std::string> tasks;
	bool isValid(const std::string & str, std::string& type, std::string& name);

public:
	OptTaskSequence(): GurlsOption(TaskSequenceOption) {};
	OptTaskSequence(const char* str): GurlsOption(TaskSequenceOption){
		tasks.push_back(str);
	}

	OptTaskSequence(std::string& str): GurlsOption(TaskSequenceOption){
		tasks.push_back(str);
	}

	OptTaskSequence& operator=(const OptTaskSequence& other);

	virtual ~OptTaskSequence(){ tasks.clear();}

	void addTask(const std::string newtask) {tasks.push_back(newtask);}
	std::vector<std::string>& getValue() {
		return tasks;
	}

	virtual bool isA(int id) { return (id == TaskSequenceOption); }
	static OptTaskSequence* dynacast(GurlsOption* opt) {
		if (opt->isA(TaskSequenceOption) ){
			return static_cast<OptTaskSequence*>(opt);
		} else {
			throw gException(gurls::Exception_Illegal_Dynamic_Cast);
		}
	}
	void getTaskAt(int index, std::string& taskdesc, std::string& taskname) {
		if (!isValid(tasks[index], taskdesc, taskname)){
			throw new gException(gurls::Exception_Invalid_TaskSequence);
		}
	}
	long int size(){
		return this->tasks.size();
	}

	//friend std::ostream& operator<<(std::ostream& os, OptString& opt);
	virtual std::ostream& operator<<(std::ostream& os);

	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int /* file_version */) const{
		ar & this->type;
		int n = this->tasks.size();
		ar & n;
		for (int i = 0; i < n; i++){
			ar & this->tasks.at(i);
		}
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int /* file_version */){
		ar & this->type;
		int n = 0;
		ar & n;
		tasks.clear();
		for(int i = 0; i < n; i++){
			std::string s("");
			ar & s;
			this->tasks.push_back(s);
		}
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()


};

/*
template <typename T>
std::ostream& operator<<(std::ostream& os, OptMatrix<T>& opt){
	return os << opt.getValue();
}*/

//template <typename T>
//std::ostream& OptMatrix<T>::operator <<(std::ostream& os){
//	return os << std::endl << this->getValue();
//}

template <typename T>
std::ostream& OptMatrix<T>::operator << (std::ostream& os){
	return os << std::endl << this->getValue();
}

//template<class Archive>
//void GurlsOption::save(Archive & ar, const unsigned int v) const {

//	int type = this->getType();
//	//ar & type;

//	switch (type) {

//	case StringOption:
//		OptString::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	case StringListOption:
//		OptStringList::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	case NumberOption:
//		OptNumber::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	case NumberListOption:
//		OptNumberList::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	case FunctionOption:
//		OptFunction::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	case MatrixOption:
//		OptMatrix<gMat2D<float> >::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	case TaskSequenceOption:
//		OptTaskSequence::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	case OptListOption:
//		GurlsOptionsList::dynacast(const_cast<GurlsOption*>(this))->save(ar, v);
//		break;
//	default:
//		break;
//	}

//}
//template<class Archive>
//void GurlsOption::load(Archive & ar, const unsigned int v){
//	int type = this->getType();
//	//ar & type;

//	switch (type)
//	{

//	case StringOption:
//		OptString::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	case StringListOption:
//		OptStringList::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	case NumberOption:
//		OptNumber::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	case NumberListOption:
//		OptNumberList::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	case FunctionOption:
//		OptFunction::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	case MatrixOption:
//		OptMatrix<gMat2D<float> >::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	case TaskSequenceOption:
//		OptTaskSequence::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	case OptListOption:
//		GurlsOptionsList::dynacast(const_cast<GurlsOption*>(this))->load(ar, v);
//		break;
//	default:
//		break;
//	}
//}



}

#endif // _GURLS_OPT_H_
