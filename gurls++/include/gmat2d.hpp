
#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4244)
#endif

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#ifdef _WIN32
#pragma warning(pop)
#endif

#include <boost/serialization/base_object.hpp>


#ifdef  USE_BINARY_ARCHIVES
typedef boost::archive::binary_iarchive iarchive;
typedef boost::archive::binary_oarchive oarchive;
#else
typedef boost::archive::text_iarchive iarchive;
typedef boost::archive::text_oarchive oarchive;
#endif

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <fstream>
#include <string>

namespace gurls {

template < typename T>
const std::string gMat2D<T>::Less = "<";
template < typename T>
const std::string gMat2D<T>::Greater = ">";
template < typename T>
const std::string gMat2D<T>::LessEq = "<=";
template < typename T>
const std::string gMat2D<T>::GreaterEq = ">=";
template < typename T>
const std::string gMat2D<T>::Equal = "==";


// IMPLEMENTATION OF TEMPLATE FUNCTIONS

template <typename T>
gMat2D<T>::gMat2D(unsigned long r, unsigned long c) : numcols(c), numrows(r)
{
    this->alloc(r*c);
}

// WARNING - HERE WE MAY BE USING A TROUBLESOME VERSION OF 'SET'
template <typename T>
gMat2D<T>::gMat2D(T* buf, unsigned long r, unsigned long c,  bool owner) : numcols(c), numrows(r)
{
    if(owner)
    {
        this->alloc(r*c);
        this->set(buf, r*c);
    }
    else
    {
        this->isowner = owner;
        this->data = buf;
    }
}

template <typename T>
gMat2D<T>::gMat2D(const gMat2D<T>& other) : numcols(other.numcols), numrows(other.numrows)
{
    this->alloc(other.numrows*other.numcols);
    *this = other;
}

// WARNING: TO BE DISCUSSED
template <typename T>
gMat2D<T>& gMat2D<T>::operator=(const gMat2D<T>& other)
{
    if (this != &other)
        this->set(other.data, other.numrows*other.numcols);

    return *this;
}

template <typename T>
void gMat2D<T>::submatrix(const gMat2D<T>& other, unsigned long r, unsigned long c)
{
    const unsigned long cols = std::min(other.numcols, this->numcols-c);
    const unsigned long rows = std::min(other.numrows, this->numrows-r);

    for(unsigned long i=0; i< cols; ++i)
        copy(this.data + r +((c+i)*this->numrows), other.data + (i*other.numrows), rows);
}

template <typename T>
void gMat2D<T>::transpose(gurls::gMat2D<T>& transposed) const
{
    if(this->numrows != transposed.numcols || this->numcols != transposed.numrows)
        throw gException(gurls::Exception_Inconsistent_Size);

    gurls::transpose(this->data, this->numrows, this->numcols, transposed.data);
}

template <typename T>
gMat2D<T> gMat2D<T>::zeros(unsigned long r, unsigned long c)
{
    gMat2D<T> m(r, c);

    gurls::set(m.data, (T)0.0, m.getSize());

    return m;
}

template <typename T>
gMat2D<T> gMat2D<T>::rand(unsigned long r, unsigned long c)
{
    gMat2D<T> m(r, c);
    T* d = m.data;

    for (unsigned long i = 0; i < r*c; ++i)
        *d++ = (T)std::rand()/RAND_MAX;

    return m;
}

template <typename T>
gMat2D<T> gMat2D<T>::eye(unsigned long n)
{
    gMat2D<T> m = gMat2D<T>::zeros(n,n);

    gurls::set(m.data, (T)1.0, n, n+1);

    return m;
}

template <typename T>
gMat2D<T> gMat2D<T>::diag(gVec<T>& diagonal)
{
    unsigned long s = diagonal.getSize();
    gMat2D<T> m = gMat2D<T>::zeros(s,s);

    gurls::copy(m.data, diagonal.getData(), s, m.numrows+1, 1);

    return m;
}


template <typename T>
void gMat2D<T>::resize(unsigned long r, unsigned long c )
{
    BaseArray<T>::resize(r*c);
    this->numrows = r;
    this->numcols = c;
}

template <typename T>
void gMat2D<T>::reshape(unsigned long r, unsigned long c){

    if (r*c != this->size)
        throw gException(gurls::Exception_Invalid_Reshape_Arguments);

    this->numcols = c;
    this->numrows = r;
}

template <typename T>
gMat2D<T> gMat2D<T>::operator+(T val) const {
    gMat2D<T> w(*this);
    w += val;
    return w;
}

/**
  * Returns a matrix containing the sum between a matrix and a scalar
  */
template <typename T>
gMat2D<T> operator+(T val, const gMat2D<T>& v) {
    gMat2D<T> w(v);
    w += val;
    return w;
}

/**
  * Returns a matrix containing the difference between a matrix and a scalar
  */
template <typename T>
gMat2D<T> operator-(T val, const gMat2D<T>& v) {
    gMat2D<T> w(-v);
    w += val;
    return w;
}

// ------------ MULT SCALAR -------------------------------------------

template <typename T>
gMat2D<T> gMat2D<T>::operator*(T val) const {
    gMat2D<T> w(*this);
    w *= val;
    return w;
}

/**
  * Returns a matrix containing the multiplication of a matrix by a scalar
  */
template <typename T>
gMat2D<T> operator*(T val, const gMat2D<T>& v) {
    gMat2D<T> w(v);
    w *= val;
    return w;
}

/**
  * Returns a matrix containing the division of a matrix by a scalar
  */
template <typename T>
gMat2D<T> operator/(T val, const gMat2D<T>& v) {
    gMat2D<T> w(v);
    w *= (static_cast<T>(1)/val);
    return w;
}

// ----------------- SUM OF VECTORS --------------------------

template <typename T>
gMat2D<T>& gMat2D<T>::operator+=(const gMat2D<T>& v) {
    return static_cast<gMat2D<T>&>( BaseArray<T>::add(v) );
}

template <typename T>
gMat2D<T> gMat2D<T>::operator+(const gMat2D<T>& v) const {
    gMat2D<T> w(v);
    w += *this;
    return w;
}


// ----------------- SUBTRACT OF VECTORS --------------------------

template <typename T>
gMat2D<T>& gMat2D<T>::operator-=(const gMat2D<T>& v) {
    return static_cast<gMat2D<T>&>( BaseArray<T>::subtract(v) );
}

template <typename T>
gMat2D<T> gMat2D<T>::operator-(const gMat2D<T>& v) const {
    gMat2D<T> w(v);
    w -= *this;
    return w;
}

// ------------------- MULT TWO VECTORS ------------------------------------

template <typename T>
gMat2D<T>& gMat2D<T>::operator*=(const gMat2D<T>& v) {
    return static_cast<gMat2D<T>&>( BaseArray<T>::multiply(v) );
}

template <typename T>
gMat2D<T> gMat2D<T>::operator*(const gMat2D<T>& v) const {
    gMat2D<T> w(v);
    w *= *this;
    return w;
}


// ------------------- DIVIDE TWO VECTORS ------------------------------------

template <typename T>
gMat2D<T>& gMat2D<T>::operator/=(const gMat2D<T>& v) {
    return static_cast<gMat2D<T>&>( BaseArray<T>::divide(v) );
}

template <typename T>
gMat2D<T> gMat2D<T>::operator/(const gMat2D<T>& v) const {
    gMat2D<T> w(v);
    w /= *this;
    return w;
}


template <typename T>
gMat2D<bool>& gMat2D<T>::operator ==(T threshold) const {
    gMat2D<bool> *w(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size; ++i, ++bptr, ++ptr){
        *bptr = (*ptr == threshold);
    }
    return *w;

}


template <typename T>
gMat2D<bool>& gMat2D<T>::operator <(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
        *bptr = (*ptr < threshold);
    }
    return *w;
}

template <typename T>
gMat2D<bool>& gMat2D<T>::operator >(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
        *bptr = (*ptr > threshold);
    }
    return *w;
}

template <typename T>
gMat2D<bool>&  gMat2D<T>::operator <=(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
        *bptr = (*ptr <= threshold);
    }
    return *w;
}

template <typename T>
gMat2D<bool>& gMat2D<T>::operator >=(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
        *bptr = (*ptr >= threshold);
    }
    return *w;
}


template <typename T>
gMat2D<bool>& gMat2D<T>::compare(T& threshold, std::string logical_operator) const {
    gMat2D<bool>* w = new gMat2D<bool>(this->rows(), this->cols());
    *w = 0;
    const T* ptr = this->data;
    bool* bptr = w->getData();

    if ( !logical_operator.compare(GreaterEq)){
        for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
            *bptr = *ptr >= threshold;
        }
    } else if (!logical_operator.compare(Greater)){
        for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
            *bptr = *ptr > threshold;
        }
    } else if (!logical_operator.compare(LessEq)){
        for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
            *bptr = *ptr <= threshold;
        }
    } else if (!logical_operator.compare(Less)){
        for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
            *bptr = *ptr < threshold;
        }
    } else if (!logical_operator.compare(Equal)){
        for(unsigned long int i=0;i<this->size;++i, ++bptr, ++ptr){
            *bptr = *ptr == threshold;
        }
    }else {
        throw gException(gurls::Exception_Logical_Operator);
    }
    return *w;
}


template <typename T>
gVec<T>& gMat2D<T>::where(const gMat2D<bool>& comparison) const {

    // WARNING: add a check on the size of the input boolean matrix

    unsigned long n0 = static_cast<unsigned long>(comparison.sum()), n1 = 0;
    //T* buf = new T[this->size];
    T* buf = new T[n0];
    const T* ptr = this->data;
    const bool* bptr = comparison.getData();
    for(unsigned long int i=0;i<this->size;++i, ++ptr, ++bptr){
        if (*bptr){
            //*(buf+i) = *ptr;
            *(buf+n1) = *ptr;
            n1++;
        }
    }
    gVec<T> *v = new gVec<T>(buf, n1, true);
    return	*v;
}



// ------ TEST IF ELEMENTS ARE ALL EQUAL TO A CONSTANT VALUE -------

template <typename T>
bool gMat2D<T>::allEqualsTo(const T& val)  const {
    const T *ptr = this->data;
    const unsigned long size = this->numrows* this->numcols;
    for (unsigned long i = 0; i < size; ++i, ++ptr){
        if (*ptr != val)
            return false;
    }
    return true;
}

// ----------------------- PRINTOUT --------------------------------

/**
  * Writes matrix information and data to a stream
  */
template <typename T>
std::ostream& operator<<(std::ostream& os, const gMat2D<T>& v)
{
   if (v.rows() >= (unsigned long)gurls::MAX_PRINTABLE_SIZE || v.cols() >= (unsigned long) gurls::MAX_PRINTABLE_SIZE)
   {
       os << v.what() << std::endl;
       return os;
   }

    os << "[";
    for (unsigned long i = 0; i < v.numrows; ++i)
    {
        for (unsigned long j = 0; j < v.numcols; ++j)
            os << " " << v.data[i+v.numrows*j];

        if( i != (v.numrows-1) )
            os << std::endl << " ";
    }

    os << " ]";
    os << std::endl;
    return os;
}


template <typename T>
void gMat2D<T>::uppertriangular(gMat2D<T>& up) const
{
    if(up.numrows != this->numrows || up.numcols != this->numcols)
        throw gException(Exception_Inconsistent_Size);

    T* upptr;
    const T* ptr;

    for (unsigned long int i = 0; i < up.numcols; ++i)
    {
        upptr = up.data + (up.numrows*(up.numcols-i-1));
        ptr = this->data + this->numrows*(this->numcols-i-1);

        gurls::set(upptr, (T)0.0, i);

        gurls::copy(upptr+i, ptr+i, up.numrows-i);
    }
}

template <typename T>
void gMat2D<T>::lowertriangular(gMat2D<T>& lo) const
{
    if(lo.numrows != this->numrows || lo.numcols != this->numcols)
        throw gException(Exception_Inconsistent_Size);

    T* loptr;
    const T* ptr;

    for (unsigned long int i = 0; i < lo.numcols; ++i)
    {
        loptr = lo.data + (lo.numrows*(lo.numcols-i-1));
        ptr = this->data + this->numrows*(this->numcols-i-1);

        gurls::copy(loptr, ptr, lo.numrows-i);

        gurls::set(loptr+lo.numrows-i, (T)0.0, i);
    }
}


template <typename T>
gMat2D<T>& gMat2D<T>::reciprocal() const {
    gMat2D<T> *w = new gMat2D(*this);
    w->setReciprocal();
    return *w;
}

template <typename T>
gVec<T>* gMat2D<T>::sum(int order) const
{
    gVec<T> *v;
    unsigned long n;

    switch(order)
    {
    case gurls::COLUMNWISE:

        n = this->numcols;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)(i).sum();

        break;
    case gurls::ROWWISE:

        n = this->numrows;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)[i].sum();

        break;
    default:
        throw gException(gurls::Exception_Illegal_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return v;
}

template <typename T>
gVec<T>* gMat2D<T>::max(int order)
{
    gVec<T> *v;
    unsigned long n;

    switch(order)
    {
    case gurls::COLUMNWISE:

        n = this->numcols;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)(i).max();

        break;
    case gurls::ROWWISE:

        n = this->numrows;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)[i].max();

        break;
    default:
        throw gException(gurls::Exception_Illegal_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return v;
}

template <typename T>
gVec<T>* gMat2D<T>::min(int order)
{
    gVec<T> *v;
    unsigned long n;

    switch(order)
    {
    case gurls::COLUMNWISE:

        n = this->numcols;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)(i).min();

        break;
    case gurls::ROWWISE:

        n = this->numrows;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)[i].min();

        break;
    default:
        throw gException(gurls::Exception_Illegal_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return v;
}

template <typename T>
gVec<T>* gMat2D<T>::argmax(int order)
{
    gVec<T> *v;
    unsigned long n;

    switch(order)
    {
    case gurls::COLUMNWISE:

        n = this->numcols;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)(i).argmax();

        break;
    case gurls::ROWWISE:

        n = this->numrows;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)[i].argmax();

        break;
    default:
        throw gException(gurls::Exception_Illegal_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return v;
}

template <typename T>
gVec<T>* gMat2D<T>::argmin(int order)
{
    gVec<T> *v;
    unsigned long n;

    switch(order)
    {
    case gurls::COLUMNWISE:

        n = this->numcols;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)(i).argmin();

        break;
    case gurls::ROWWISE:

        n = this->numrows;
        v = new gVec<T>(n);

        for (unsigned long i = 0; i < n; ++i)
            v->at(i) = (*this)[i].argmin();

        break;
    default:
        throw gException(gurls::Exception_Illegal_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return v;
}

template <typename T>
template<class Archive>
void gMat2D<T>::save(Archive & ar, const unsigned int /* file_version */) const
{
    ar & this->numrows & this->numcols & this->isowner;

    T* ptr = this->data;
    T* ptr_end = this->data+(this->numrows*this->numcols);

    while (ptr!=ptr_end)
        ar & *ptr++;
}

template <typename T>
template<class Archive>
void gMat2D<T>::load(Archive & ar, const unsigned int /* file_version */)
{
    ar & this->numrows;
    ar & this->numcols;
    ar & this->isowner;

    this->size = this->numrows*this->numcols;
    this->data = new T[this->size];
    T* ptr = this->data;
    T* ptr_end = this->data+this->size;

    while (ptr!=ptr_end)
        ar & *ptr++;
}

template <typename T>
void gMat2D<T>::load(const std::string& fileName)
{
#ifndef USE_BINARY_ARCHIVES
    std::ifstream instream(fileName.c_str());
#else
    std::ifstream instream(fileName.c_str(), std::ios_base::binary);
#endif

    if(!instream.is_open())
        throw gException("Could not open file " + fileName);

    try
    {
        iarchive inar(instream);
        inar >> *this;
    }
    catch(boost::archive::archive_exception&)
    {
        instream.close();
        throw gException("Invalid file format for " + fileName);
    }

    instream.close();
}

template <typename T>
void gMat2D<T>::save(const std::string& fileName) const
{
#ifndef USE_BINARY_ARCHIVES
    std::ofstream outstream(fileName.c_str());
#else
    std::ofstream outstream(fileName.c_str(), std::ios_base::binary);
#endif

    if(!outstream.is_open())
        throw gException("Could not open file " + fileName);

    oarchive outar(outstream);
    outar << *this;

    outstream.close();
}

template <typename T>
void gMat2D<T>::readCSV(const std::string& fileName, bool colMajor)
{
    std::vector<std::vector< T > > matrix;
    std::ifstream in(fileName.c_str());

    if(!in.is_open())
        throw gurls::gException("Cannot open file " + fileName);

    unsigned long rows = 0;
    unsigned long cols = 0;

    std::string line;
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(" ;|,");

    while (std::getline(in, line))
    {
        if(!line.empty())
        {
            std::vector<T> tf;

            tokenizer tokens(line, sep);
            for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it)
                tf.push_back(boost::lexical_cast<T>(*it));

            matrix.push_back(tf);
            ++rows;
        }
    }

    in.close();

    if(!matrix.empty())
        cols = matrix[0].size();

    this->resize(rows, cols);

    if(colMajor)
    {
        for(unsigned long i=0; i<rows; ++i)
            gurls::copy(this->data + i, &(*(matrix[i].begin())), cols, rows, 1);
    }
    else
    {
        for(unsigned long i=0; i<rows; ++i)
            gurls::copy(this->data + i*cols, &(*(matrix[i].begin())), rows);
    }
}

template <typename T>
void gMat2D<T>::saveCSV(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());

    if(!out.is_open())
        throw gurls::gException("Cannot open file " + fileName);

    for (unsigned long i = 0; i < numrows; ++i)
    {
        for (unsigned long j = 0; j < numcols; ++j)
            out << " " << this->data[i+numrows*j];

        out << std::endl;
    }

    out.close();
}

}
