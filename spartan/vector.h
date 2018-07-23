// -*- C++-*-

#ifndef __SPARTAN_vector_h__
#define __SPARTAN_vector_h__

#include <iostream>
#include <cmath>

namespace spartan {

  class Vector {
  public:    
    Vector(int size = 0) ;
    Vector(int size , double t) ;
    ~Vector() ;
    
    double& operator[](int index) ;
    const double& operator[](int index) const ;
    
    inline int size() const {return size_;}
    inline void resize(int n) {allocate(n);} 
    void normalize() ;
    void zero() ;

    inline double* data(){return data_;}
    inline const double* data() const {return data_;}
    
    Vector& operator=(const Vector& A) ;    
    Vector& operator+=(const Vector& A)  ;
    Vector& operator-=(const Vector& A)  ;
    Vector& operator*=(const double& scal)  ; 

  private:    
    double* data_ ;
    int size_ ;
    void allocate(int size) {
      if(data_) delete[] data_ ; 
      data_ = (size == 0) ? NULL : new double[size]; 
      size_ = size ;
    }      
  };
  
  inline double norm2( const Vector& x) {
    double result = 0.0; 
    for(int i=0; i<x.size(); ++i){
      result += std::pow(x[i],2) ;
    }
    return result ;
  }  

  inline double norm( const Vector& x) {
    return std::sqrt(norm2(x)); 
  }  

  inline double dot( const Vector& x, 
	      const Vector& y){
    double result = 0.0; 
    for(int i=0; i<x.size(); ++i){
      result += x[i]*y[i] ;
    }
    return result ;
  } 

  // y = y + a*x
  inline void axpy(Vector& y,  
	    double a, 
	    const Vector& x ) {
    for(int i=0; i<x.size(); ++i){
      y[i] += a*x[i] ;
    }
  } 
  
  // x = a*x
  inline void scal(Vector& x, 
	    const double& a ){
    for(int i=0; i<x.size(); ++i){
      x[i] *= a ; 
    }
  } 

  // x = a*x +y
  inline void scal_add(Vector& x, 
		const double& a,
		const Vector& y){
    for(int i=0; i<x.size(); ++i){
      x[i] = y[i] + a*x[i] ; 
    }
  } 
  
  // y = x
  inline void copy( const Vector& x, Vector& y ) {
    for(int i=0; i<x.size(); ++i){
      y[i] = x[i] ; 
    }
  } 
  /*    
  
  */
}; 

inline std::ostream& operator<<(std::ostream& o, const spartan::Vector& x) {
    int n = x.size(); 
    if(n>0){
      o<< "("<<n<<") [";
      for(int i=0; i<n-1; ++i){
	o<< x[i] <<" " ;
      }
      o<<x[n-1]<<"]"<<std::endl ; 

    }
    return o ; 
  }

#endif
