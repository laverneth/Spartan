#include <iostream>
#include "spartan/vector.h"
#include "spartan/sparse_types.h"

using namespace spartan ;


int main(){
  // 1D Laplacian
  int n = 1000 ;
  Matrix A ;
  A.n = n ;
  A.nnz = 3*n-2;
  A.v = new double[3*n-2] ; 
  A.irow = new int[n+1] ;
  A.jcol = new int[3*n-2];
  int nnz = 0 ; 
  Vector x(n);
  for(int i=0; i<n; ++i){
    A.irow[i] = nnz ; 
    x[i] = 1.0;
    // Sub-diagonal
    if(i>0){
      A.v[nnz] = -1.0 ;
      A.jcol[nnz] = i-1;
      nnz++; 
    }
    // diagonal
    A.v[nnz] = 2.0 ;
    A.jcol[nnz] = i ; 
    nnz++;
    // Sup-diagonal
    if(i+1<n){
      A.v[nnz] = -1.0 ;
      A.jcol[nnz] = i+1 ;
      nnz++; 
    } 
  }
  A.irow[n] = nnz ; 

  std::cout<<" Matrix is build !" << std::endl ; 

  std::cout<<" x " << x <<std::endl ; 
  
}
