#include "gurls++/traintest.h"

namespace gurls
{
//"easy train" fucntion
//requires in input:
//		- X, a T(float or double) buffer, intended as a matrix n x d
//		- Y, a T(float or double) buffer, intended as a matrix n x t
//		- n, d and t, geometry parameters described above
//		- algorithm, an optional string containing the selected algorithm (only "krls" available now)
//		- kernel, an optional string containing the selected kernel type ("gaussian" and "linear" available now)
//		- problem, an optional string containing the selected problem type ("classification" and "regression" available)
//		- savefile, an optional string containing the path in which to save the model
//gives in output
//		- a GurlsOptionList with the trained model, ready to be used with "test" function
template <typename T>
gurls::GurlsOptionsList train(T* X, T* y, unsigned long n, unsigned long d, unsigned long t, 
							std::string algorithm, std::string kernel, std::string problem, std::string savefile)
{
    gMat2D<T> Xtr(X, n, d, false);
    gMat2D<T> ytr(y, n, t, false);

    try
    {
		GurlsWrapper<T> *gurlsWrap;
		
		if(algorithm=="krls")
			gurlsWrap= new KernelRLSWrapper<T>("Train");
		else{
			algorithm="krls";
			gurlsWrap= new KernelRLSWrapper<T>("Train");
			std::cout<<"algorithm forced to krls"<<std::endl;}
		
		if(kernel=="linear" && algorithm=="krls")
			 dynamic_cast<KernelRLSWrapper<T>*>(gurlsWrap)->setKernelType(KernelWrapper<T>::LINEAR);
		else if(algorithm=="krls")
			dynamic_cast<KernelRLSWrapper<T>*>(gurlsWrap)->setKernelType(KernelWrapper<T>::RBF);

		if(savefile!="")
			gurlsWrap->setSavefile(savefile);	
		else
			gurlsWrap->setSavefile("temp");

		if(problem=="classification"){
			gurlsWrap->setProblemType(gurlsWrap->CLASSIFICATION);}
		else if(problem=="regression")
			gurlsWrap->setProblemType(gurlsWrap->REGRESSION);
		else{
			const char * probList[] = { "Classification", "Regression"};
            typename KernelWrapper<T>::ProblemType prob =gurlsWrap->problemTypeFromData(Xtr, ytr);
			std::cout<<"Problem set to "<<probList[prob];
			gurlsWrap->setProblemType(prob);}

		gurlsWrap->train(Xtr, ytr);
		if(savefile=="")
			std::remove("temp");
		GurlsOptionsList retopt=gurlsWrap->getOpt();
		
		delete gurlsWrap;
		return retopt;
    }

    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        throw e;
    }
}
//"easy test" fucntion
//requires in input:
//		- model, previously saved opt, already trained on data using "train" function
//		- X, a T(float or double) buffer, intended as a matrix n x d
//		- predbuff, a T(float or double) buffer, intended as a matrix n x t, where prediction will be saved
//		- Y, an optional T(float or double) buffer, intended as a matrix n x t
//		- perfbuff, an optional T(float or double) buffer, intended as a t length vector, where performance will be saved
//		- n, d and t, geometry parameters described above
//		- perfstring, a string containing the type of perf to be calculated, "rmse", "macroavg" and "auto" are the accepted values,
//					  leave blank to skip performance calculation (perfbuff will remain untouched)
//gives in output
//		- EXIT_SUCCESS or EXIT_FAILURE
template <typename T>
int test(gurls::GurlsOptionsList& model, T* X, T* Y, T* predbuff, T* perfbuff, unsigned long n, unsigned long d, unsigned long t, std::string perfstring)
{
	
    try
    {
	GurlsWrapper<T> *gurlsWrap;
	if (model.hasOpt("kernel"))
		{
			gurlsWrap = new KernelRLSWrapper<T>("test");
			std::string kernelType = model.getOptValue<OptString>("kernel.type");
			if (kernelType=="linear")
				dynamic_cast<KernelRLSWrapper<T>*>(gurlsWrap)->setKernelType(KernelWrapper<T>::LINEAR);
		}
	else
		{
			std::cout<<"no kernel found... setting to linear";
			gurlsWrap = new KernelRLSWrapper<T>("test");
			dynamic_cast<KernelRLSWrapper<T>*>(gurlsWrap)->setKernelType(KernelWrapper<T>::LINEAR);
		}

	gurlsWrap->loadOpt(model);
	gMat2D<T> *perfMat;
    gMat2D<T> Xte(X, n, d, false);

	if(perfstring!="macroavg" && perfstring != "rmse" && perfstring != "" && t>0 || perfstring=="auto")
	{
		if(model.hasOpt("hoperf"))
			perfstring=model.getOptAsString("hoperf");
		else
			perfstring="macroavg";
		std::cout<<"Performance measure set to default for selected problem: "<<perfstring<<std::endl;
	}

        //load the test data
	gMat2D<T> *pred= gurlsWrap->eval(Xte);
	copy(predbuff, pred->getData(), pred->getSize());

	if (t>0 && perfstring != "")
	{
		gMat2D<T> yte(Y, n, t, false);
		perfMat = gurlsWrap->perf(yte, *pred, perfstring);
		copy(perfbuff, perfMat->getData(), perfMat->getSize());
		delete perfMat;
	}
	
	delete gurlsWrap;
    return EXIT_SUCCESS;
    }
    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}
//"easy test" function, with load from file
//requires in input:
//		- loadfile, path to previously saved optionlist in .bin format, already trained on data using "train" function
//		- X, a T(float or double) buffer, intended as a matrix n x d
//		- predbuff, a T(float or double) buffer, intended as a matrix n x t, where prediction will be saved
//		- Y, an optional T(float or double) buffer, intended as a matrix n x t
//		- perfbuff, an optional T(float or double) buffer, intended as a t length vector, where performance will be saved
//		- n, d and t, geometry parameters described above
//		- perfstring, a string containing the type of perf to be calculated, "rmse", "macroavg" and "auto" are the accepted values,
//					  leave blank to skip performance calculation (perfbuff will remain untouched)
//gives in output
//		- EXIT_SUCCESS or EXIT_FAILURE
template <typename T>
int test(std::string loadfile, T* X, T* Y, T* predbuff, T* perfbuff, unsigned long n, unsigned long d, unsigned long t, std::string perfstring)
{
	
	GurlsOptionsList model("model");
	model.load(loadfile);

    try
    {
	GurlsWrapper<T> *gurlsWrap;
	if (model.hasOpt("kernel"))
		{
			gurlsWrap = new KernelRLSWrapper<T>("test");
			std::string kernelType = model.getOptValue<OptString>("kernel.type");
			if (kernelType=="linear")
				dynamic_cast<KernelRLSWrapper<T>*>(gurlsWrap)->setKernelType(KernelWrapper<T>::LINEAR);
		}

	gurlsWrap->loadOpt(model);
	gMat2D<T> *perfMat;
    gMat2D<T> Xte(X, n, d, false);

	if(perfstring!="macroavg" && perfstring != "rmse" && perfstring != "" && t>0 || perfstring=="auto")
	{
		if(model.hasOpt("hoperf"))
			perfstring=model.getOptAsString("hoperf");
		else
			perfstring="macroavg";
		std::cout<<"Performance measure set to default for selected problem: "<<perfstring<<std::endl;
	}

        //load the test data
	gMat2D<T> *pred= gurlsWrap->eval(Xte);
	copy(predbuff, pred->getData(), pred->getSize());

	if (t>0 && perfstring != "")
	{
		gMat2D<T> yte(Y, n, t, false);
		perfMat = gurlsWrap->perf(yte, *pred, perfstring);
		copy(perfbuff, perfMat->getData(), perfMat->getSize());
		delete perfMat;
	}
	
	delete gurlsWrap;
    return EXIT_SUCCESS;
    }
    catch (gException& e)
    {
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}

}
