
namespace gurls {

// IMPLEMENTATION OF TEMPLATE METHODS

template <typename T>
BaseArray<T>::BaseArray(const BaseArray<T>& other) {
    // TO BE DISCUSSED: do we need to allocate memory here?
    // this->alloc(other.size);
    *this = other;
}

template <typename T>
BaseArray<T>& BaseArray<T>::operator=(const BaseArray<T>& other) {
    if (this == &other) {
        return *this;
    } else {
        this->alloc(other.size);
        this->set(other.data, other.size);
    }
}

template <typename T>

BaseArray<T>& BaseArray<T>::operator=(const T& val) {
    T* end = &(this->data[this->size]);
    for (T* ptr = this->data; ptr != end; ++ptr) {
        *ptr = val;
    }
    return *this;
}


template <typename T>
void BaseArray<T>::alloc(unsigned long n) {
    this->isowner = true;
    this->size = n;
    if ( this->size > 0 )
        this->data = new T[this->size];
}

template <typename T>
void BaseArray<T>::set(const T* buf, unsigned long n, unsigned long start) {
    // An exception should be raised here instead of using assert
    assert(n <= (this->size - start));
    std::memcpy((this->data)+start, buf, n*sizeof(T));
}


template <typename T>
void BaseArray<T>::randomize(){
    T* d = this->data;
    for (unsigned long i = 0; i < this->size; i++) {
        *d++ = (T)std::rand()/RAND_MAX;
    }
}

template <typename T>
void BaseArray<T>::resize(unsigned long n) {
    if (! this->isowner) {
        throw gException("It is not possible to resize an array that is not the owner of its data.");
    }
    else if (n != this->size) {
        T* tmp = this->data;
        unsigned long oldsize = this->size;
        this->alloc(n);
        this->set(tmp, std::min(n, oldsize));
        if(tmp != NULL)
            delete[] tmp;
    };
}

template <typename T>
void BaseArray<T>::asarray(T * buf, unsigned long n) const {
    if (n > this->size) {
        throw gException("The length of the array is smaller than the required number of elements to be copied.");
    }
    std::memcpy(buf, this->data, n*sizeof(T));
}


/*
  template <typename T>
  BaseArray<T> BaseArray<T>::operator-() const {
  //return static_cast<T>(-1)*(*this);
  BaseArray<T> copy(*this);
  return copy *= -1;
  }
*/

template <typename T>
BaseArray<T>& BaseArray<T>::operator+=(T val) {
    T* ptr = this->data;
    for (unsigned long i = 0; i < this->size; ++i, ++ptr) {
        *ptr += val;
    }
    return *this;
}

template <typename T>
BaseArray<T>& BaseArray<T>::operator-=(T val) {
    return *this += (-val);
}

template <typename T>
BaseArray<T>& BaseArray<T>::operator*=(T val) {
    T* ptr = this->data;
    for (unsigned long i = 0; i < this->size; ++i, ++ptr) {
        *ptr *= val;
    }
    return *this;
}

template <typename T>
BaseArray<T>& BaseArray<T>::operator/=(T val) {
    T* ptr = this->data;
    for (int i = 0; i < this->size; ++i, ++ptr) {
        *ptr /= val;
    }
    return *this;
}

template <typename T>
BaseArray<T>& BaseArray<T>::add(const BaseArray<T>& v) {
    T *ptr = this->data, *ptr_v = v.data;
    for (unsigned long i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
        *ptr += *ptr_v;
    }
    return *this;
}

template <typename T>
BaseArray<T>& BaseArray<T>::subtract(const BaseArray<T>& v) {
    T *ptr = this->data, *ptr_v = v.data;
    for (int i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
        *ptr -= *ptr_v;
    }
    return *this;
}

template <typename T>
BaseArray<T>& BaseArray<T>::multiply(const BaseArray<T>& v) {
    T *ptr = this->data, *ptr_v = v.data;
    for (unsigned long i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
        *ptr *= *ptr_v;
    }
    return *this;
}

template <typename T>
BaseArray<T>& BaseArray<T>::divide(const BaseArray<T>& v) {
    T *ptr = this->data, *ptr_v = v.data;
    for (unsigned long i = 0; i < this->size; ++i, ++ptr, ++ptr_v) {
        *ptr /= *ptr_v;
    }
    return *this;
}

/**
  * Checks if all elements in a vector are equal to a given value
  */
template <typename U>
bool operator== (const BaseArray<U>& v, const U& val) {
    U *ptr = v.data;
    for (unsigned int i = 0; i < v.size; ++i, ++ptr) {
        if (*ptr != val){
            return false;
        }
    }
    return true;
}

template <typename T>
bool BaseArray<T>::closeTo(const BaseArray<T>& other, T tolerance) const {
    if (this->size != other.size) {
        return false;
    }
    T *ptr0 = this->data, *ptr1 = other.data;
    for (unsigned int i = 0; i < this->size; ++ptr0, ++ptr1, ++i) {
        if (std::fabs(*ptr0 - *ptr1) > tolerance) {
            return false;
        }
    }
    return true;
}

template <typename T>
const T& BaseArray<T>::max() const {
    if (this->size == 0){return *new T(std::numeric_limits<T>::max());}
    const T* val = std::max_element(&(this->data[0]), &(this->data[this->size]));
    return *val;
}

template <typename T>
const T& BaseArray<T>::min() const {
    if (this->size == 0){return *new T(std::numeric_limits<T>::min());}
    const T* val = std::min_element(&(this->data[0]), &(this->data[this->size]));
    return *val;
}

template <typename T>
double BaseArray<T>::sum() const{
    T* ptr = this->data;
    T* endptr = this->data+this->size;

    double sum = 0.0;
    while(ptr != endptr){
        sum+=*ptr++;
    }
    return sum;
}


template <typename T>
BaseArray<T>& BaseArray<T>::setReciprocal() {
    T* ptr = this->data;
    T* endptr = this->data+this->size;
    T one = static_cast<T>(1.0);
    while(ptr != endptr){
        *ptr = one / *ptr;
        ++ptr;
    }
    return *this;
}

template <typename T>
template<class Archive>
void BaseArray<T>::save(Archive & ar, const unsigned int /* file_version */) const{
    ar & this->size & this->isowner;
    T* ptr = this->data;
    T* ptr_end = this->data+this->size;
    while (ptr!=ptr_end){
        ar & *ptr++;
    }
}

template <typename T>
template<class Archive>
void BaseArray<T>::load(Archive & ar, const unsigned int /* file_version */){
    ar & this->size;
    ar & this->isowner;
    this->data = new T[this->size];
    T* ptr = this->data;
    T* ptr_end = this->data+this->size;
    while (ptr!=ptr_end){
        ar & *ptr++;
    }
}

}
