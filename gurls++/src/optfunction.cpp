#include "optfunction.h"

#include <ostream>

namespace gurls
{


OptFunction::OptFunction(): GurlsOption(FunctionOption), f(NULL){}


OptFunction::OptFunction(std::string func_name): GurlsOption(FunctionOption)
{
    setValue(func_name);
}

OptFunction::~OptFunction()
{
    if(f != NULL)
        delete f;
}

OptFunction& OptFunction::operator=(const OptFunction& other)
{
    if(f != NULL)
        delete f;

    setValue(other.name);
    return *this;
}

std::string OptFunction::getName() const
{
    return name;
}

bool OptFunction::isA(OptTypes id) const
{
    return (id == FunctionOption);
}

OptFunction* OptFunction::dynacast(GurlsOption* opt)
{
    if (opt->isA(FunctionOption) )
        return static_cast<OptFunction*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

const OptFunction* OptFunction::dynacast(const GurlsOption* opt)
{
    if (opt->isA(FunctionOption) )
        return static_cast<const OptFunction*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

GURLS_EXPORT std::ostream& OptFunction::operator<<(std::ostream& os)
{
    os << "Pointer to the function <" << this->getName()
       << "> whose signature is: T (*func)(T*, int)" ;
    return os;
}

void OptFunction::setValue(std::string func_name)
{
    name = func_name;

    if(func_name == "mean")
        f = new Mean();
    else if(func_name == "min")
        f = new Min();
    else if(func_name == "max")
        f = new Max();
    else if(func_name == "median")
        f = new Median();
    else
    {
        f = NULL;
        throw gException(Exception_Unknown_Function);
    }
}

}
