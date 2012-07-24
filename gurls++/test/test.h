#ifndef TEST_H
#define TEST_H

#include <exception>
#include <string>
#include <list>
#include <utility>
#include <fstream>

#include "gmat2d.h"
#include "optlist.h"
#include "exceptions.h"
#include "gmath.h"

#include "rlsprimalr.h"
#include "rlsdualr.h"

#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

//#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace boost::filesystem;

namespace gurls
{

namespace test
{

template<typename T>
void check_vector(const T* result, const T* reference, const unsigned long size)
{
    const T errCoeff = 1.0e-10;

//    std::cout.precision(8);

//    for(const T *res_it = result, *res_end = res_it+size, *ref_it = reference; res_it != res_end; ++res_it, ++ref_it)
//        std::cout << (*res_it) << "   " << (*ref_it) << std::endl;

    for(const T *res_it = result, *res_end = res_it+size, *ref_it = reference; res_it != res_end; ++res_it, ++ref_it)
        BOOST_REQUIRE_LE(std::abs(*res_it - *ref_it), errCoeff*std::abs(*ref_it));
}

template<>
void check_vector(const unsigned long* result, const unsigned long* reference, const unsigned long size)
{
//    for(const unsigned long *res_it = result, *res_end = res_it+size, *ref_it = reference; res_it != res_end; ++res_it, ++ref_it)
//        std::cout << (*res_it) << "   " << (*ref_it) << std::endl;

    for(const unsigned long *res_it = result, *res_end = res_it+size, *ref_it = reference; res_it != res_end; ++res_it, ++ref_it)
        BOOST_REQUIRE_EQUAL(*res_it - *ref_it, *ref_it);
}

template<typename T>
void check_matrix(gurls::gMat2D<T>& result, gurls::gMat2D<T>& reference)
{
//    std::cout << result.rows() << " " << reference.rows() << std::endl;
//    std::cout << result.cols() << " " << reference.cols() << std::endl;

    BOOST_REQUIRE_EQUAL(result.rows(), reference.rows());
    BOOST_REQUIRE_EQUAL(result.cols(), reference.cols());

    check_vector(result.getData(), reference.getData(), result.getSize());
}

template<typename T>
gMat2D<T> * readMatrix(const std::string &fileName, bool ROWM = true )
{
    std::vector<std::vector< double > > matrix;
    std::ifstream in(fileName.c_str());

    if(!in.is_open())
        throw gurls::gException("Cannot open file " + fileName);

    unsigned long rows = 0;
    unsigned long cols = 0;

    std::string line;
    while (std::getline(in, line))
    {
        std::istringstream ss(line);
        std::vector<double> tf;
        std::copy(std::istream_iterator<double>(ss), std::istream_iterator<double>(), std::back_inserter(tf));

        matrix.push_back(tf);
        ++rows;
    }
    in.close();

    if(matrix.empty())
        cols = 0;
    else
        cols = matrix[0].size();

    gMat2D<T> *ret =  new gMat2D<T>(rows, cols);
    T* buffer = ret->getData();

    for(int i=0; i<rows; ++i)
    {
        for(int j=0; j<cols; ++j)
        {
            if(ROWM)
                buffer[i*cols+j]= static_cast<T>(matrix[i][j]);
            else
                buffer[j*rows+i]= static_cast<T>(matrix[i][j]);
        }
    }

    return ret;
}

template<typename T>
GurlsOption* openFile(std::string fileName, OptTypes type)
{
    std::ifstream file(fileName.c_str());

    switch(type)
    {
    case StringOption:
    {
        file.seekg (0, std::ios::end);
        int length = file.tellg();
        file.seekg (0, std::ios::beg);

        char* data = new char[length+1];

        file.read(data, length);
        data[length] = '\0';

        file.close();

        OptString* ret = new OptString(data);
        delete[] data;
        return ret;
    }
    case NumberOption:
    {
        double ret = 0.0;

        file >> ret;
        file.close();

        return new OptNumber(ret);
    }

    case StringListOption:
    {
        std::vector<std::string> ret;

        std::copy(std::istream_iterator<std::string>(file), std::istream_iterator<std::string>(), std::back_inserter(ret));

        file.close();
        return new OptStringList(ret);
    }
    case NumberListOption:
    {
        std::vector<double> ret;

        std::copy(std::istream_iterator<double>(file), std::istream_iterator<double>(), std::back_inserter(ret));

        file.close();
        return new OptNumberList(ret);
    }
    case FunctionOption:
    {
        std::string ret;

        file >> ret;

        file.close();
        return new OptFunction(ret);
    }

    case MatrixOption:
    case VectorOption:
        file.close();
        return new OptMatrix<gMat2D<T> >(*(readMatrix<T>(fileName)));

    case GenericOption:
    case OptListOption:
    case TaskSequenceOption:
    case TaskIDOption:
        file.close();
        throw "Unsupported option";
    }
}

template<typename T>
void checkOptions(GurlsOption& result, GurlsOption& reference)
{
    BOOST_REQUIRE_EQUAL( result.getType(), reference.getType());

    switch(result.getType())
    {
    case StringOption:
        BOOST_REQUIRE_EQUAL (OptString::dynacast(&result)->getValue(), OptString::dynacast(&reference)->getValue());
        break;
    case NumberOption:
        check_vector<double>( &(OptNumber::dynacast(&result)->getValue()), &(OptNumber::dynacast(&reference)->getValue()), 1);
        break;
    case StringListOption:
    {
        std::vector<std::string>& res = OptStringList::dynacast(&result)->getValue();
        std::vector<std::string>& ref = OptStringList::dynacast(&reference)->getValue();

        BOOST_REQUIRE_EQUAL(res.size(), ref.size());

        for(int i=0; i< res.size(); ++i)
            BOOST_REQUIRE_EQUAL(res[i], ref[i]);

        break;
    }
    case NumberListOption:
    {
        std::vector<double>& res = OptNumberList::dynacast(&result)->getValue();
        std::vector<double>& ref = OptNumberList::dynacast(&reference)->getValue();

        BOOST_REQUIRE_EQUAL(res.size(), ref.size());

        check_vector<double>( &(*(res.begin())), &(*(ref.begin())), res.size());

        break;
    }
    case FunctionOption:
        BOOST_REQUIRE_EQUAL(OptFunction::dynacast(&result)->getName(), OptFunction::dynacast(&reference)->getName());
        break;

    case MatrixOption:
    case VectorOption:
        check_matrix(OptMatrix<gMat2D<T> >::dynacast(&result)->getValue(), OptMatrix<gMat2D<T> >::dynacast(&reference)->getValue());
        break;

    case GenericOption:
    case OptListOption:
    case TaskSequenceOption:
    case TaskIDOption:
        throw "Unsupported option";
    }
}

class Data
{
public:
    Data(path dataDirectory, path taskName, bool loadTrainingData)
    {
        path dataPath = dataDirectory / taskName;

        if(loadTrainingData)
        {
            X = (dataDirectory / path("Xtr.txt")).native();
            Y = (dataDirectory / path("ytr.txt")).native();
        }
        else
        {
            X = (dataDirectory / path("Xte.txt")).native();
            Y = (dataDirectory / path("yte.txt")).native();
        }

        readIndexFile(dataPath / path("in.txt"), in);
        readIndexFile(dataPath / path("out.txt"), out);
    }

    void loadDefaults()
    {
        defaults.push_back(PairType("hoperf.txt", StringOption));
        defaults.push_back(PairType("singlelambda.txt", FunctionOption));
        defaults.push_back(PairType("nlambda.txt", NumberOption));
        defaults.push_back(PairType("nsigma.txt", NumberOption));
        defaults.push_back(PairType("nholdouts.txt", NumberOption));
        defaults.push_back(PairType("hoproportion.txt", NumberOption));
        defaults.push_back(PairType("smallnumber.txt", NumberOption));
    }

    typedef std::pair<std::string, OptTypes> PairType;
    typedef std::list<PairType> ListType;

    std::string X;
    std::string Y;
    ListType defaults;
    ListType in;
    ListType out;

protected:
    OptTypes string2OptTypes(std::string type)
    {
        if(type == "matrix")
            return MatrixOption;
        if(type == "string")
            return StringOption;
        if(type == "number")
            return NumberOption;
        if(type == "function")
            return FunctionOption;
        if(type == "stringlist")
            return StringListOption;
        if(type == "numberlist")
            return NumberListOption;

        assert(0);
        return GenericOption;
    }

    void readIndexFile(path filePath, ListType& list)
    {
        if(!exists(filePath))
            throw gurls::gException(std::string("Cannot open file ") + filePath.native());

        std::ifstream stream( filePath.native().c_str());

        std::string name;
        std::string type;

        while(stream.good())
        {
            stream >> name >> type;

            if(stream.good())
                list.push_back(PairType(name, string2OptTypes(type)));
        }
        stream.close();
    }

};


template<typename T, class TaskType>
class Fixture
{
public:
    Fixture(path dataDirectory, path taskName, Data& _data) : dataDir(dataDirectory), task(taskName), data(_data)
    {
        typedef typename Data::ListType ListType;
        typedef typename Data::PairType PairType;

        path dataPath = this->dataDir / this->task;

        X = readMatrix<T>(data.X);
        Y = readMatrix<T>(data.Y);

        opt = new gurls::GurlsOptionsList("Testdata");

        for(typename ListType::iterator it = data.defaults.begin(), end = data.defaults.end(); it != end; ++it)
        {
            path fileName(it->first);
            path filePath(path(dataDir / fileName));

            if(!exists(filePath))
                throw gurls::gException(std::string("Cannot open file ") + filePath.native());

            GurlsOption* option = openFile<T>(filePath.native(), it->second);
            opt->addOpt(fileName.stem().native(), option);
        }

        for(typename ListType::iterator it = this->data.in.begin(), end = this->data.in.end(); it != end; ++it)
        {
            PairType& p = *it;
            path fileName(p.first);
            path filePath(path(dataPath / fileName));

            if(!exists(filePath))
                throw gurls::gException(std::string("Cannot open file ") + filePath.native());

            std::vector<std::string> strs;
            boost::split(strs, fileName.stem().native(), boost::is_any_of("-"));

            GurlsOption* option;

            if(fileName.stem() == "split-indices" || fileName.stem() == "split-lasts")
                option = openFile<unsigned long>(filePath.native(), MatrixOption);
            else if(fileName.stem() == "paramsel-lambdas")
                option = openFile<T>(filePath.native(), NumberListOption);
            else
                option = openFile<T>(filePath.native(), p.second);

            GurlsOptionsList* opt_it = this->opt;
            for(int i=0; i< strs.size()-1; ++i)
            {
                if( !(opt_it->hasOpt(strs[i])) )
                {
                    GurlsOptionsList* tmp = new GurlsOptionsList(strs[i]);
                    opt_it->addOpt(strs[i], tmp);
                    opt_it = tmp;
                }
                else
                    opt_it = GurlsOptionsList::dynacast(opt_it->getOpt(strs[i]));
            }

            opt_it->addOpt(strs.back(), option);
        }
    }

    ~Fixture()
    {
        delete X;
        delete Y;
        delete opt;
    }

    void runTask()
    {
        TaskType taskProcess;
        taskProcess.execute(*X, *Y, *opt);
    }

    void checkResults(std::string optField)
    {
        typedef typename Data::ListType ListType;
        typedef typename Data::PairType PairType;

        path dataPath = this->dataDir / this->task;

        GurlsOption* fieldToCheck = this->opt->getOpt(optField);

        for(typename ListType::iterator it = this->data.out.begin(), end = this->data.out.end(); it != end; ++it)
        {
            PairType& p = *it;
            path fileName(p.first);
            path filePath(path(dataPath / fileName));

            if(!exists(filePath))
                throw gurls::gException(std::string("Cannot open file ") + filePath.native());

            if(fileName.stem() == "perf-forplot")
                continue;

            std::vector<std::string> strs;
            boost::split(strs, fileName.stem().native(), boost::is_any_of("-"));

            BOOST_REQUIRE_EQUAL( strs[0], optField);

            GurlsOption* result;

            if(strs.size() > 1)
            {
                GurlsOptionsList* opt_it = GurlsOptionsList::dynacast(fieldToCheck);

                for(int i=1; i< strs.size()-1; ++i)
                    opt_it = GurlsOptionsList::dynacast(opt_it->getOpt(strs[i]));

                result = opt_it->getOpt(strs.back());
            }
            else
                result = fieldToCheck;

            GurlsOption* reference;

//            std::cout << fileName.stem() << std::endl << std::endl;

            if(fileName.stem() == "split-indices" || fileName.stem() == "split-lasts")
            {
                reference = openFile<unsigned long>(filePath.native(), MatrixOption);
                checkOptions<unsigned long>(*result, *reference);
            }
            else
            {
                if(fileName.stem() == "paramsel-lambdas")
                    reference = openFile<T>(filePath.native(), NumberListOption);
                else
                    reference = openFile<T>(filePath.native(), p.second);

                checkOptions<T>(*result, *reference);
            }

            delete reference;
        }
    }

    gurls::gMat2D<T>* X;
    gurls::gMat2D<T>* Y;
    gurls::GurlsOptionsList* opt;

protected:

    path dataDir;
    path task;

    Data& data;

};

}

}

#endif //TEST_H
