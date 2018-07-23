#include "preconditioner.h"

namespace spartan{

  

  DiagonalPreconditioner::DiagonalPreconditioner(const Matrix& A){
    diag_.resize(A.n);
    for(int i=0; i<A.n; ++i){
      for(int j=A.irow[i]; j<A.irow[i+1]; ++j){
	if(A.jcol[j] == i){
	  diag_[i] = 1.0/A.v[j] ; 
	}
      }
    }
  }


  void DiagonalPreconditioner::apply(const Vector& x, Vector& y) const {
    for(int i=0; i<x.size(); ++i){
      y[i] = x[i]*diag_[i] ; 
    }
  }
  
}
