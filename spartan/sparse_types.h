// -*- C++ -*-
//

#ifndef __SPARTAN_sparse_types_h__
#define __SPARTAN_sparse_types_h__

#include "vector.h"

namespace spartan {
  
   struct Coeff {
      double value ;
      int index ;
   } ;

  class SparseArray {
  public:
    SparseArray() ;
    ~SparseArray() ;
    
    int size() const ;

    Coeff & coeff(int ii) ;
    const Coeff& coeff(int ii) const ;

    void add(int index, double val) ;
    void sort() ;
  
    void clear() ;
    void zero() ;

    bool value(int ii, double& value) const;
 
  protected:
    void grow() ;

  private:    
    Coeff* coeff_ ;
    // size is the actual size of values which are stored
    int size_ ;
    // capacity is the maximum number of values which can be stored
    int capacity_ ;
  } ;

  class TempMatrix {
    

  }; 
  
  /**
   * @brief standard class of matrix, stored in CRS
   */
  struct Matrix {
    Matrix(); 
    ~Matrix();
    void reset(); 
    int n;
    int nnz ; 
    int* irow ;
    int* jcol; 
    double* v ;
    bool is_symmetric;

    void mult(const Vector& x,
	      Vector& y) const ; 
    void mult_transpose(const Vector& x,
			Vector& y) const ;

  };

  

  
  struct Cholesky {
    Matrix ld;
    Matrix ldT;
    int* p;
    int* invp; 
  };
  
  static void cuthill_mac_kee(Matrix& A, int* perm) ;


  inline int SparseArray::size() const { 
    return size_ ; 
  }  
  
  inline Coeff& SparseArray::coeff(int ii) { 
    return coeff_[ii] ;
  }  
  
  inline const Coeff& SparseArray::coeff(int ii) const {
    return coeff_[ii] ;
  }  
  
}
#endif
