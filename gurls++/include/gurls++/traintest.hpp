#include "gurls++/traintest.h"

namespace gurls
{

template <typename T>
gurls::GurlsOptionsList train(T* X, T* y, unsigned long n, unsigned long d, unsigned long t, 
							std::string algorithm, std::string kernel, std::string problem, std::string savefile)
{
    gMat2D<T> Xtr(X, n, d, false);
    gMat2D<T> ytr(y, n, t, false);
	GurlsWrapper<T> *gurlsWrap;

    try
    {
		if(algorithm=="krls")
			gurlsWrap= new KernelRLSWrapper<T>("Train");
		else{
			algorithm="krls";
			gurlsWrap= new KernelRLSWrapper<T>("Train");
			std::cout<<"algorithm forced to krls"<<std::endl;
			}
		
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
			std::cout<<"Problem automatically set to: "<<probList[prob];
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
		delete gurlsWrap;
        throw e;
    }
}

template <typename T>
int test(gurls::GurlsOptionsList& model, T* X, T* Y, T* predbuff, T* perfbuff, unsigned long n, unsigned long d, unsigned long t, std::string perfstring)
{
	
	GurlsWrapper<T> *gurlsWrap;
    try
    {
	if (model.hasOpt("kernel"))
		{
			gurlsWrap = new KernelRLSWrapper<T>("test");
			std::string kernelType = model.getOptValue<OptString>("kernel.type");
			if (kernelType=="linear")
				dynamic_cast<KernelRLSWrapper<T>*>(gurlsWrap)->setKernelType(KernelWrapper<T>::LINEAR);
		}
	else
		{
			std::cout<<"No kernel found... setting to linear"<<std::endl;
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
	}
	std::cout<<"Performance measure set to: "<<perfstring<<std::endl;

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
		delete gurlsWrap;
        return EXIT_FAILURE;
    }
}

template <typename T>
int test(std::string loadfile, T* X, T* Y, T* predbuff, T* perfbuff, unsigned long n, unsigned long d, unsigned long t, std::string perfstring)
{
	
	GurlsOptionsList model("model");
	model.load(loadfile);
	GurlsWrapper<T> *gurlsWrap;

    try
    {
	if (model.hasOpt("kernel"))
		{
			gurlsWrap = new KernelRLSWrapper<T>("test");
			std::string kernelType = model.getOptValue<OptString>("kernel.type");
			if (kernelType=="linear")
				dynamic_cast<KernelRLSWrapper<T>*>(gurlsWrap)->setKernelType(KernelWrapper<T>::LINEAR);
		}
	else
		{
			std::cout<<"No kernel found... setting to linear"<<std::endl;
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
	}
	std::cout<<"Performance measure set to: "<<perfstring<<std::endl;
	
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
		delete gurlsWrap;
        std::cout << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}

}
