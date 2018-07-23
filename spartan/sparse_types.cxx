#include "sparse_types.h"
#include "allocators.h"
#include  <algorithm>

namespace spartan {

  SparseArray::SparseArray() {
    coeff_ = Allocator<Coeff>::allocate(2) ;
    size_ = 0 ;
    capacity_ = 2 ;
  }
    
  SparseArray::~SparseArray() {
    Allocator<Coeff>::deallocate(coeff_) ;
  }
  
 
  inline void SparseArray::clear() { 
    Allocator<Coeff>::deallocate(coeff_) ;
    coeff_ = Allocator<Coeff>::allocate(2) ;
    size_ = 0 ;
    capacity_ = 2 ;
  }  
  
  inline void SparseArray::zero() {
    size_ = 0 ;
  }
  
  void SparseArray::add(int index, double val) {
    Coeff* coeff = NULL ;
    // search if coefficient already exist
    for(int ii=0; ii < size_; ii++) {
      if(coeff_[ii].index == index) {
	coeff = &(coeff_[ii]) ;
	break ;
      }
    }
    
    if(coeff != NULL) {
      coeff->value += val ;
    } 
    else {
      size_++ ;
      if(size_ > capacity_) {
	grow() ;
      }
      coeff = &(coeff_[size_ - 1]) ;
      coeff->value = val ;
      coeff->index = index ;
    }
  }
  
  void SparseArray::grow() {
    int old_capacity = capacity_ ;
    if(old_capacity<8){
      capacity_ = capacity_ *2 ;
    }
    else {
      capacity_ += 8 ; 
    }
    Allocator<Coeff>::reallocate(coeff_, old_capacity, capacity_) ;
  }
  
  
  class CoeffIndexCompare {
  public:
    bool operator()(const Coeff& c1, 
		    const Coeff& c2) {
      return c1.index < c2.index ;
    }
  } ;
  
  
  void SparseArray::sort() {
    Coeff* begin = coeff_ ;
    Coeff* end   = coeff_ + size_ ;
    std::sort(begin, end, CoeffIndexCompare()) ;
  }  
  
  bool SparseArray::value(int index, double& value) const {
    for(int ii=0; ii < size_; ii++) {
      const Coeff& coeff = coeff_[ii];
      if(coeff.index == index) {
	value = coeff.value ;
	return true ;
      }
    }
    value = 0 ;
    return false ;
  }
  


  Matrix::Matrix(){
    // constructor of arrays is not systematically called
    reset() ; 
  }

  Matrix::~Matrix(){
    // destructor of arrays is not systematically called
    reset() ; 
  }

  void Matrix::reset(){
    n = -1 ;
    nnz = -1; 
    irow  = nullptr; 
    jcol = nullptr; 
    v = nullptr;
  }

  void Matrix::mult(const Vector& x,
		    Vector& y) const {
    
    for(int i=0; i<n; ++i){
      double val = 0.0; 
      for(int j= irow[i]; j< irow[i+1]; ++j) {
	val += v[j]*x[jcol[j]] ; 
      }
      y[i] = val ; 
    }

  }
  
  void Matrix::mult_transpose(const Vector& x,
			      Vector& y) const {

 
    for(int i=0; i<n; ++i) y[i] = 0.0; 
    for(int i=0; i<n; ++i){
      const double& xi = x[i];
      for(int j=irow[i]; j<irow[i+1]; ++j)
	y[jcol[j]] += v[j]*xi;
    } 
  }


  void cuthill_mac_kee(Matrix& A, int* perm) {
    // step 1: find a node with lowest degree
    int lowest_degree_index = -1;
    {
      int lowest_degree = A.n;
      for(int i=0; i<A.n; ++i){
	for(int j= A.irow[i]; j< A.irow[i+1]; ++j) {
	  int nj = A.irow[i+1] - A.irow[i] ; 
	  if(nj<lowest_degree){
	    lowest_degree = nj ;
	    lowest_degree_index = j;
	  }
	}
      }
    }
    // step 2 bfs traversal of matrix graph
    
  }
}
