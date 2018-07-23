#include "spartan/vector.h"
#include "spartan/sparse_types.h"
#include "spartan/pcg.h"
#include "spartan/preconditioner.h"
#include "spartan/hsc_preconditioner.h"
#include <fstream>
#include <vector>
#include <map>

using namespace spartan ;


static void load(const std::string& filename, Matrix& crs) {	
  std::cout<<" load file: "<<filename<< std::endl ;
  std::ifstream myfile (filename);
  std::vector<std::map<int, double> > Arep ; 
  int nrow, ncol, nnz ;
  int sym_nnz  = 0 ; 
  if (myfile.is_open()){    
    myfile>>nrow>>ncol>>nnz ;
    Arep.resize(nrow); 
    std::cout<<" nrow: "<< nrow << std::endl ;
    std::cout<<" nnz: "<< nnz << std::endl ;
    int i, j ;
    double val;
    for(int k=0; k<nnz; ++k){
      myfile>>i>>j>>val ;
      i--;
      j--;
      Arep[i][j] = val ;
      sym_nnz++;
      if(j != i){
	Arep[j][i] = val ; 
	sym_nnz++;
      }
    }
    myfile.close();
  }
  std::cout<<"nnz "<<nnz<<" sym: "<<sym_nnz<< std::endl ; 
  crs.n = nrow ;
  crs.v = new double[sym_nnz] ;
  crs.irow = new int[nrow+1];
  crs.jcol = new int[sym_nnz] ; 

  // convert Arep to crs
  nnz = 0;
  for(int i=0; i<nrow; ++i){
    crs.irow[i]= nnz ;
    for(std::map<int, double>::iterator it = Arep[i].begin(); it != Arep[i].end(); ++it){
      crs.v[nnz]  = it->second;
      crs.jcol[nnz] = it->first ;
      ++nnz ; 
    }
  }
  crs.irow[nrow]= nnz ;
}
  



int main(){

  Matrix A ; 
  load("/home/tom/Desktop/SparseMatrix/Spartan/data/msc01050.mtx", A);
  int n = A.n ;
  HSCPreconditioner prec(A); 
  
  
  /*
  Vector x(n, 0.0);
  Vector b(n, 1.0);
  int nb_pcg_iter = pcg(A, prec, b, x, 5000, 1.e-16); 
  std::cout<<" nb iter: "<< nb_pcg_iter << std::endl ; 
  Vector res(n, 0.0);
  A.mult(x, res);
  axpy(res, -1, b); 
  std::cout<<" residual norm: "<< norm2(res) << std::endl;
  */
}; 
