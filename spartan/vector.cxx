#include "vector.h"


namespace spartan {
   
  Vector::Vector(int size ) {

    data_ = NULL ;
    allocate(size) ;
  } 

  Vector::Vector(int size , double t) {
    data_ = NULL ;
    allocate(size) ;
    for(int i=0; i<size; ++i) data_[i] = t ;
  } 

  Vector::~Vector() {
    delete[] data_ ;
    data_ = NULL ;
  }
    
  double& Vector::operator[](int index) {
    return data_[index] ;
  }
    
  const double& Vector::operator[](int index) const {
    return data_[index] ;
  }
  
  void Vector::normalize() {
    double norm = std::sqrt(norm2(*this)); 
    if(norm>1.e-20){
      for(int i=0; i<size_; ++i)
	data_[i] /= norm ;
    }
  }

  void Vector::zero() {
    for(int i=0; i<size_; ++i) data_[i] = double(0) ; 
  }

  Vector& Vector::operator=(const Vector& A) {
    if (size_ == A.size_){
      // do nothing !
    }
    else{
      allocate(A.size_);
    }
    for(int i=0; i<size_; ++i){
      data_[i] = A.data_[i];
    }      
    return *this;
  }
    
  Vector& Vector::operator+=(const Vector& A)  {
    for(int i=0; i<size_; i++){
      data_[i] += A.data_[i];
    }
    return (*this);
  }

  Vector& Vector::operator-=(const Vector& A)  {
    for(int i=0; i<size_; i++){
      data_[i] -= A.data_[i];
    }
    return (*this);
  }

  Vector& Vector::operator*=(const double& scal)  {
    for(int i=0; i<size_; i++){
      data_[i] *= scal;
    }
    return (*this);
  }
        
}
