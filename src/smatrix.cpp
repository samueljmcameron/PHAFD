#include <complex>
#include <stdexcept>

#include "smatrix.hpp"

using namespace PHAFD_NS;

template <typename T>
sMatrix<T>& sMatrix<T>::operator= (const sMatrix & rhs)
{
  if (&rhs != this) {

    uncreate();

    create(rhs.begin(), rhs.end());

    _sizeax.at(0) = rhs.axis_size(0);
    _sizeax.at(1) = rhs.axis_size(1);
  }

  return *this;

}

template <typename T>
void sMatrix<T>::create()
{
  arr = limit = 0;
}


template <typename T>
void sMatrix<T>::create(size_type n, const T& val)
{
  arr = alloc.allocate(n);
  limit = arr + n;
  std::uninitialized_fill(arr, limit,val);

}

template <typename T>
void sMatrix<T>::create(const_iterator i, const_iterator j)
{

  arr = alloc.allocate(j - i);
  limit = std::uninitialized_copy(i,j,arr);
  
}

template <typename T>
void sMatrix<T>::uncreate()
{
  if (arr) {

    iterator it = limit;
    while (it != arr)
      alloc.destroy(--it);

    alloc.deallocate(arr,limit-arr);
  }

  arr = limit = 0;
}


template <typename T>
void sMatrix<T>::resize(size_type Ny, size_type Nx)
{

  size_type new_size = Ny*Nx;

  iterator new_arr = alloc.allocate(new_size);
  iterator new_limit = std::uninitialized_copy(arr, limit, new_arr);

  uncreate();

  arr = new_arr;
  limit = new_limit;

  _sizeax.at(0) = Ny;
  _sizeax.at(1) = Nx;
  
  
}

template class sMatrix<double>;
template class sMatrix<std::complex<double>>;
