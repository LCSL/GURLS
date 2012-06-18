
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/base_object.hpp>


#ifdef  USE_BINARY_ARCHIVES
typedef boost::archive::binary_iarchive iarchive;
typedef boost::archive::binary_oarchive oarchive;
#else
typedef boost::archive::text_iarchive iarchive;
typedef boost::archive::text_oarchive oarchive;
#endif

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
gMat2D<T>::gMat2D(unsigned long r, unsigned long c) : numcols(c), numrows(r) {
    this->alloc(r*c);
}

// WARNING - HERE WE MAY BE USING A TROUBLESOME VERSION OF 'SET'
template <typename T>
gMat2D<T>::gMat2D(T* buf, unsigned long r, unsigned long c,  bool owner) : numcols(c), numrows(r) {
    if (owner) {
        this->alloc(r*c);
        this->set(buf, r*c);
    } else {
        this->isowner = owner;
        this->data = buf;
    };
}

template <typename T>
gMat2D<T>::gMat2D(const gMat2D<T>& other) : numcols(other.numcols), numrows(other.numrows) {
    this->alloc(other.numrows*other.numcols);
    *this = other;
}

// WARNING: TO BE DISCUSSED
template <typename T>
gMat2D<T>& gMat2D<T>::operator=(const gMat2D<T>& other) {
    if (this != &other) {
        this->set(other.data, other.numrows*other.numcols);
    }
    return *this;
}

template <typename T>
void gMat2D<T>::submatrix(const gMat2D<T>& other, unsigned long r, unsigned long c) {
    unsigned long w1 = this->cols();
    unsigned long w2 = other.cols();
    unsigned long size = std::min(w2, w1-c);
    const T* ptrin = other.getData();
    int i1, i2;
    for (int i = 0; i < other.rows() && r+i < this->rows(); ++i) {
        i1 = (i+r)*w1;
        i2 = i*w2;
        for (int j = 0; j < w2 && c+j < w1; ++j) {
            this->data[i1+j+c] = ptrin[i2+j];
        }
    }
}

template <typename T>
void gMat2D<T>::transpose(gMat2D <T>& transposed) const {
    if (this->rows() != transposed.cols() || this->cols() != transposed.rows()) {
        throw gException(gurls::Exception_Inconsistent_Size);
    }
    T* d1 = transposed.getData();
    T* d0 = this->data;
    int N = this->cols();

    for (unsigned long c = 0; c < this->cols(); ++c) {
        for (unsigned long r = 0; r < this->rows(); ++r){
            *d1++=*(d0+r*N+c);
        }
    }
}


template <typename T>
gMat2D<T> gMat2D<T>::zeros(unsigned long r, unsigned long c) {
    gMat2D<T> m(r, c);
    memset(m.data, 0, r*c*sizeof(T));
    return m;
}

template <typename T>
gMat2D<T> gMat2D<T>::rand(unsigned long r, unsigned long c) {
    gMat2D<T> m(r, c);
    T* d = m.getData();
    for (unsigned long i = 0; i < r*c; i++) {
        *d++ = (T)std::rand()/RAND_MAX;
    }
    return m;
}

template <typename T>
gMat2D<T> gMat2D<T>::eye(unsigned long n){
    gMat2D<T> m = gMat2D<T>::zeros(n,n);
    T* d = m.getData();
    T one = static_cast<T>(1.0);
    for (int i = 0; i < n; i++) {
        d[i*(n+1)] = one;
    }
    return m;
}


template <typename T>
gMat2D<T> gMat2D<T>::diag(gVec<T> diagonal) {
    unsigned long s = diagonal.getSize();
    gMat2D<T> m = gMat2D<T>::zeros(s,s);
    T* d = m.getData();
    for (int i = 0; i < s; i++) {
        d[i*(s+1)] = diagonal[i];
    }
    return m;
}


// WARNING: TO BE DISCUSSED
template <typename T>
void gMat2D<T>::resize(unsigned long r, unsigned long c ) {
    BaseArray<T>::resize(r*c);
    this->numrows = r;
    this->numcols = c;
}

template <typename T>
void gMat2D<T>::reshape(unsigned long r, unsigned long c){

    if (r*c != this->size) {
        throw gException(gurls::Exception_Invalid_Reshape_Arguments);
    }
    this->numcols = c;
    this->numrows = r;

}



template <typename T>
gMat2D<T> gMat2D<T>::operator+(T val) const {
    gMat2D<T> w(*this);
    w += val;
    return w;
}

template <typename T>
gMat2D<T> operator+(T val, const gMat2D<T>& v) {
    gMat2D<T> w(v);
    w += val;
    return w;
}

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

template <typename T>
gMat2D<T> operator*(T val, const gMat2D<T>& v) {
    gMat2D<T> w(v);
    w *= val;
    return w;
}


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

    for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
        *bptr = (*ptr == threshold);
    }
    return *w;

}


template <typename T>
gMat2D<bool>& gMat2D<T>::operator <(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
        *bptr = (*ptr < threshold);
    }
    return *w;
}

template <typename T>
gMat2D<bool>& gMat2D<T>::operator >(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
        *bptr = (*ptr > threshold);
    }
    return *w;
}

template <typename T>
gMat2D<bool>&  gMat2D<T>::operator <=(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
        *bptr = (*ptr <= threshold);
    }
    return *w;
}

template <typename T>
gMat2D<bool>& gMat2D<T>::operator >=(T threshold) const {
    gMat2D<bool> *w = new gMat2D<bool>(this->rows(), this->cols());
    const T* ptr = this->data;
    bool* bptr = w->getData();

    for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
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
        for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
            *bptr = *ptr >= threshold;
        }
    } else if (!logical_operator.compare(Greater)){
        for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
            *bptr = *ptr > threshold;
        }
    } else if (!logical_operator.compare(LessEq)){
        for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
            *bptr = *ptr <= threshold;
        }
    } else if (!logical_operator.compare(Less)){
        for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
            *bptr = *ptr < threshold;
        }
    } else if (!logical_operator.compare(Equal)){
        for(unsigned long int i=0;i<this->size;i++, bptr++, ptr++){
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
    for(unsigned long int i=0;i<this->size;i++, ptr++, bptr++){
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

template <typename T>
std::ostream& operator<<(std::ostream& os, const gMat2D<T>& v) {
   if (v.rows() >= (unsigned long)gurls::MAX_PRINTABLE_SIZE || v.cols() >= (unsigned long) gurls::MAX_PRINTABLE_SIZE){
       os << v.what() << std::endl;
       return os;
   }
    os << "[";
    for (unsigned long i = 0; i < v.numrows; ++i){
        for (unsigned long j = 0; j < v.numcols; ++j) {
            os << " " << v.data[i*v.numcols+j];
        }
        if( i == (v.numrows-1) )
            os << " ]";
        else{
            os << std::endl << " ";
        }
    }
    os << std::endl;
    return os;
}


template <typename T>
void gMat2D<T>::uppertriangular(gMat2D<T>& up) const{

    T* upptr = up.getData();
    const T* ptr = this->data;
    for (unsigned long int i = 0; i < up.rows(); i++){
        for (unsigned long int j = 0; j < i; j++, ptr++, upptr++){
            *upptr=0;
        }
        for (unsigned long int j = i; j < up.cols(); j++, ptr++, upptr++) {
            *upptr = *ptr;
        }
    }

}

template <typename T>
void gMat2D<T>::lowertriangular(gMat2D<T>& lo) const {

    T* upptr = lo.getData();
    const T* ptr = this->data;
    for (unsigned long int i = 0; i < lo.rows(); i++){
        for (unsigned long int j = 0; j <= i; j++, ptr++, upptr++){
            *upptr = *ptr;
        }
        for (unsigned long int j = i+1; j < lo.cols(); j++, ptr++, upptr++) {
            *upptr=0;
        }
    }

}


template <typename T>
gMat2D<T>& gMat2D<T>::reciprocal() const {
    gMat2D<T> *w = new gMat2D(*this);
    w->setReciprocal();
    return *w;
}

template <typename T>
gVec<T>& gMat2D<T>::sum(int order) {
    gVec<T> *v;
    unsigned long n;
    if (order == gurls::COLUMNWISE){
        n = this->cols();
        v = new gVec<T>(n);
        *v = 0;
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = static_cast<T>((this->operator ()(i)).sum());
        }
    }else if (order == gurls::ROWWISE){
        n = this->rows();
        v = new gVec<T>(n);
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = static_cast<T>((this->operator [](i)).sum());
        }
    } else {
        throw gException(gurls::Exception_Illegat_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return *v;
}

template <typename T>
gVec<T>& gMat2D<T>::max(int order) {
    gVec<T> *v;
    unsigned long n;
    if (order == gurls::COLUMNWISE){
        n = this->cols();
        v = new gVec<T>(n);
        *v = 0;
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator ()(i)).max();
        }
    }else if (order == gurls::ROWWISE){
        n = this->rows();
        v = new gVec<T>(n);
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator [](i)).max();
        }
    } else {
        throw gException(gurls::Exception_Illegat_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return *v;
}

template <typename T>
gVec<T>& gMat2D<T>::min(int order) {
    gVec<T> *v;
    unsigned long n;
    if (order == gurls::COLUMNWISE){
        n = this->cols();
        v = new gVec<T>(n);
        *v = 0;
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator ()(i)).min();
        }
    }else if (order == gurls::ROWWISE){
        n = this->rows();
        v = new gVec<T>(n);
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator [](i)).min();
        }
    } else {
        throw gException(gurls::Exception_Illegat_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return *v;
}



template <typename T>
gVec<T>& gMat2D<T>::argmax(int order) {
    gVec<T> *v;
    unsigned long n;
    if (order == gurls::COLUMNWISE){
        n = this->cols();
        v = new gVec<T>(n);
        *v = 0;
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator ()(i)).argmax();
        }
    }else if (order == gurls::ROWWISE){
        n = this->rows();
        v = new gVec<T>(n);
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator [](i)).argmax();
        }
    } else {
        throw gException(gurls::Exception_Illegat_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return *v;
}

template <typename T>
gVec<T>& gMat2D<T>::argmin(int order) {
    gVec<T> *v;
    unsigned long n;
    if (order == gurls::COLUMNWISE){
        n = this->cols();
        v = new gVec<T>(n);
        *v = 0;
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator ()(i)).argmin();
        }
    }else if (order == gurls::ROWWISE){
        n = this->rows();
        v = new gVec<T>(n);
        for (unsigned long i = 0; i < n; i++){
            v->at(i) = (this->operator [](i)).argmin();
        }
    } else {
        throw gException(gurls::Exception_Illegat_Argument_Value+" -> order must be either gurls::COLUMNWISE or gurls::ROWWISE.");
    }

    return *v;
}

template <typename T>
template<class Archive>
void gMat2D<T>::save(Archive & ar, const unsigned int /* file_version */) const{
    ar & this->numrows & this->numcols & this->isowner;
    T* ptr = this->data;
    T* ptr_end = this->data+(this->numrows*this->numcols);
    while (ptr!=ptr_end){
        ar & *ptr++;
    }
}

template <typename T>
template<class Archive>
void gMat2D<T>::load(Archive & ar, const unsigned int /* file_version */){
    ar & this->numrows;
    ar & this->numcols;
    ar & this->isowner;
    this->size = this->numrows*this->numcols;
    this->data = new T[this->size];
    T* ptr = this->data;
    T* ptr_end = this->data+this->size;
    while (ptr!=ptr_end){
        ar & *ptr++;
    }
}

template <typename T>
void gMat2D<T>::load(const std::string& fileName)
{
    std::ifstream instream(fileName.c_str());

    if(!instream.is_open())
        throw gException("Could not open file " + fileName);

    iarchive inar(instream);
    inar >> *this;

    instream.close();
}

template <typename T>
void gMat2D<T>::save(const std::string& fileName) const
{
    std::ofstream outstream(fileName.c_str());

    if(!outstream.is_open())
        throw gException("Could not open file " + fileName);

    oarchive outar(outstream);
    outar << *this;

    outstream.close();
}



}
