// -*- C++ -*-

#ifndef __SPARTAN_pcg_h__
#define __SPARTAN_pcg_h__

#include "sparse_types.h"
#include "vector.h"
#include "preconditioner.h"

namespace spartan {

  int pcg(const Matrix& A, 
	  const Preconditioner& prec,
	  const Vector& b, 
	  Vector& x, 
	  int max_iter, 
	  double eps) {
    
    int N = A.n;     
    
    if( N==0 ){
      return 0 ;
    }
    
    Vector r(N) ;
    Vector d(N) ; 
    Vector h(N) ; 
    Vector& Ad = h ; 
  
    int iter=0;
    double rh, alpha, beta;
    double err = eps*eps*dot(b, b) ; 
    double r_r ; 

    r.zero() ; 
    x.zero();
    
    // r =  Ax
    A.mult(x,r);

    // r = r - b 
    axpy(r, -1.0, b) ;
    prec.apply(r,d);

    // h = d 
    copy(d,h) ;     
    
    rh = dot(r,h);
    r_r = dot(r,r); 
    
    while ( (r_r>err) && (iter<max_iter) ){
    
      A.mult(d,Ad);
      alpha = rh/dot(d,Ad);      
      axpy(x, -alpha, d) ; 
      axpy(r, -alpha, Ad) ; 

      // apply preconditioner
      prec.apply(r,h);

      beta = 1./rh; 
      rh = dot(r,h);
      beta *= rh ;
      d *= beta ; 
      axpy(d, 1.0, h) ; 

      // recompute residual
      r_r = dot(r,r);
      ++iter; 
    }
    return iter;
  }

}

#endif
