#ifndef PSPDE_SMATRIX_HPP
#define PSPDE_SMATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <memory>

namespace PHAFD_NS {
template <typename T>
class sMatrix
{
public:

  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  
  sMatrix() : _sizeax(2)
  {
    create();
    _sizeax.at(0) = 0;
    _sizeax.at(1) = 0;
  };
  explicit sMatrix(size_type Ny, size_type Nx, const T& t = T())
    : _sizeax(2)
  {
    create(Nx*Ny,t);
    _sizeax.at(0) = Ny;
    _sizeax.at(1) = Nx;
  }

  sMatrix (const sMatrix & m)
    : _sizeax(2)
  { // copy constructor
    
    create(m.begin(), m.end());
    _sizeax.at(0) = m.axis_size(0);
    _sizeax.at(1) = m.axis_size(1);
    
  } 
  sMatrix& operator= (const sMatrix &);
  

  size_type size() const { return limit - arr; }

  
  T& operator()(size_type i, size_type j) { return arr[i*_sizeax[1] + j]; }
  const T& operator()(size_type i, size_type j) const { return arr[i*_sizeax[1] + j]; }

  iterator begin() { return arr; }
  const_iterator begin() const { return arr; }


  iterator end() { return limit; }
  const_iterator end() const { return limit; }

  ~sMatrix() { uncreate(); };
  
  
  T* data() {
      return arr;
  };

  void resize(size_type, size_type);
  
  int axis_size(int i) const
  {
    return _sizeax.at(i);
  }

  int size_including(int);

  
private:
  iterator arr;
  iterator limit;

  std::vector<size_type> _sizeax;

  std::allocator<T> alloc;
  

  void create();
  void create( size_type, const T&);
  void create(const_iterator, const_iterator);

  void uncreate();

};




};
#endif
