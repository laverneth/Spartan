// -*- C++ -*-
//

#ifndef __SPARTAN_allocator_h__
#define __SPARTAN_allocator_h__

#include "vector.h"

namespace spartan {

  /**
   * @brief contiguous array-type of allocator
   */
  template <class T>
  class Allocator {
  public:
    static inline T*   allocate(int number) ;
    static inline void deallocate(T*& addr) ;
    static inline void reallocate(T*& addr, int old_nbr, int new_nbr) ;
  } ;
  
  template <class T>
  inline T* Allocator<T>::allocate(int number) {
    return new T[number] ;
  }
  
  template <class T>
  inline void Allocator<T>::deallocate(T*& addr) {
    if(addr) delete[] addr ;
    addr = NULL ; 
  }
  
  template <class T>
  inline void Allocator<T>::reallocate(T*& addr, int old_nbr, int new_nbr) {
    T* new_addr = new T[new_nbr] ;
    for(int i=0; i<old_nbr; i++) {
      new_addr[i] = addr[i] ;
    }
    delete[] addr ;
    addr = new_addr ;
  }

}
#endif
