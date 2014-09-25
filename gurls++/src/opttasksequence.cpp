#include "gurls++/opttasksequence.h"

namespace gurls
{


bool OptTaskSequence::isValid(const std::string & str, std::string& type, std::string& name)
{
    size_t found = str.find(gurls::TASKDESC_SEPARATOR);

    if (found==std::string::npos)
        return false;

    type = str.substr(0, found);
    name = str.substr(found+1);

    if (name.find(gurls::TASKDESC_SEPARATOR)!=std::string::npos)
        return false;

    return true;
}

void OptTaskSequence::addTask(const std::string newtask)
{
    m_tasks.push_back(new OptTask(newtask));
}

bool OptTaskSequence::isA(OptTypes id) const
{
    return (id == TaskSequenceOption);
}

OptTaskSequence *OptTaskSequence::dynacast(GurlsOption *opt)
{
    if (opt->isA(TaskSequenceOption) )
        return static_cast<OptTaskSequence*>(opt);

    throw gException(gurls::Exception_Illegal_Dynamic_Cast);
}

const OptTaskSequence *OptTaskSequence::dynacast(const GurlsOption *opt)
{
    if (opt->isA(TaskSequenceOption) )
        return static_cast<const OptTaskSequence*>(opt);

    throw gException(Exception_Illegal_Dynamic_Cast);
}

unsigned long OptTaskSequence::size() const
{
    return (unsigned long)m_tasks.size();
}

void OptTaskSequence::clear()
{
    m_tasks.clear();
}


}
