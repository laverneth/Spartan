#include "spartan/vector.h"
#include "spartan/sparse_types.h"
#include "spartan/pcg.h"
#include "spartan/preconditioner.h"
#include <fstream>
#include <vector>
#include <map>

using namespace spartan ;


static void sym_load(const std::string& filename, Matrix& crs) {	
  std::cout<<" sym load file: "<<filename<< std::endl ;
  std::ifstream myfile (filename);
  int nrow, ncol, nnz ;
  if (myfile.is_open()){    
    myfile>>nrow>>ncol>>nnz ;
    crs.n = nrow ;
    crs.v = new double[nnz] ;
    crs.irow = new int[nrow+1];
    crs.jcol = new int[nnz] ;
    std::cout<<" nrow: "<< nrow << std::endl ;
    std::cout<<" nnz: "<< nnz << std::endl ;
    double kmem = (((nrow+1 +nnz)/double(1024))*sizeof(int)+(nnz/double(1024))*sizeof(double)) ;
    double Mmem = kmem/1024. ;
    std::cout<<" estimated size: "<< int(Mmem) <<" Mb " <<std::endl;
    int i, j ;
    double val;
    int ilast = -1 ; 
    int nz = 0 ;
    for(int k=0; k<nnz; ++k){
      myfile>>i>>j>>val ;
      i--;
      //if(i% (nrow/10) == 0) std::cout<<" load: "<< i/(nrow/10)<<std::endl ; 
      j--;
      if(i != ilast){
	crs.irow[i]= nz ;
	ilast = i;
      }
      crs.v[nz]  = val; 
      crs.jcol[nz] = j; 
      ++nz ; 
    }
    crs.irow[nrow]= nz ;
  }  
  myfile.close();
}


static void total_load(const std::string& filename, Matrix& crs) {	
  std::cout<<" sym load file: "<<filename<< std::endl ;
  std::ifstream myfile (filename);
  int n, ncol, nnz ;
  if (myfile.is_open()){    
    myfile>>n>>ncol>>nnz ;

    SparseArray* tmp_rows = new SparseArray[n] ;  
    std::cout<<" n  : "<< n << std::endl ;
    std::cout<<" nnz: "<< nnz << std::endl ;
    double kmem = (((n+1 +nnz)/double(1024))*sizeof(int)+ (nnz/double(1024))*sizeof(double)) ;
    double Mmem = kmem/1024. ;
    std::cout<<" estimated size: "<< int(Mmem) <<" Mb " <<std::endl;
    int i, j ;
    double val;
    
    for(int k=0; k<nnz; ++k){
      myfile>>i>>j>>val ;
      i--;
      j--;
      tmp_rows[i].add(j, val) ;
      if(j != i){
	tmp_rows[j].add(i, val) ;
      }
    }

    // sort
    for(int i=0; i<n; ++i){
      tmp_rows[i].sort() ;
    }
    // Fill the Matrix
    int sym_nnz = 2*(nnz-n) + n ; 
    crs.n = n ;
    crs.v = new double[sym_nnz] ;
    crs.irow = new int[n+1];
    crs.jcol = new int[sym_nnz] ;
    

  }
  
  myfile.close();
}

static void load(const std::string& filename, Matrix& crs) {	
  std::cout<<" load file: "<<filename<< std::endl ;
  std::ifstream myfile (filename);
  std::vector<std::map<int, double> > Arep ; 
  int nrow, ncol, nnz ;

  int sym_nnz  = 0 ; 
  if (myfile.is_open()){    
    myfile>>nrow>>ncol>>nnz ;
    std::cout<<" estimated size: "<< (nrow*sizeof(int)+nnz*sizeof(int) +nnz*sizeof(double))/(1024*1024)<<" Mb " <<std::endl;
    Arep.resize(nrow); 
    std::cout<<" nrow: "<< nrow << std::endl ;
    std::cout<<" nnz: "<< nnz << std::endl ;
    int i, j ;
    double val;
    for(int k=0; k<nnz; ++k){
      myfile>>i>>j>>val ;
      i--;
      if(i% (nrow/10) == 0) std::cout<<" load: "<< i/(nrow/10)<<std::endl ; 
      j--;
      Arep[i][j] = val ;
      sym_nnz++;
      /*if(j != i){
	Arep[j][i] = val ; 
	sym_nnz++;
      }
      */
    }
    myfile.close();
  }
  std::cout<<"Matrix is loaded in temporary format! nnz "<<nnz<<" sym: "<<sym_nnz<< std::endl ; 
  crs.n = nrow ;
  crs.v = new double[sym_nnz] ;
  crs.irow = new int[nrow+1];
  crs.jcol = new int[sym_nnz] ; 

  // convert Arep to crs
  std::cout<<" Convert matrix to symmetric CRS \n";
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
  //load("/home/tom/Desktop/SparseMatrix/Spartan/data/msc01050.mtx", A);
  //load("/home/tom/Desktop/SparseMatrix/Spartan/data/fv3.mtx", A);
  //load("/home/tom/Desktop/SparseMatrix/Spartan/data/thermomech_TC.mtx", A);
  sym_load("/home/tom/Desktop/SparseMatrix/Spartan/data/bone010.mtx", A);
  int n = A.n ;
  std::cout<<" problem size " << n <<std::endl;
  int row_size_min = n ;
  int row_size_max = 0;
  double row_size_avg = 0; 
  for(int i=0; i<n; ++i){
    int row_size_i = A.irow[i+1] - A.irow[i];
    row_size_min = std::min(row_size_min, row_size_i); 
    row_size_max = std::max(row_size_max, row_size_i); 
    row_size_avg += double(row_size_i) ; 
  }
  row_size_avg /= double(n); 
  std::cout<<" row size min: "<< row_size_min << std::endl ;
  std::cout<<" row size max: "<< row_size_max << std::endl ;
  std::cout<<" row size avg: "<< row_size_avg << std::endl ;
  /*
  std::cout<<" Compute N-score "<<std::endl;
  double nscore = 0.0 ;
  for(int i=0; i<n; ++i){
    int jmax = -1;
    for(int j = A.irow[i]; j<A.irow[i+1]; ++j){
      if(A.jcol[j]>jmax)
	jmax = A.jcol[j];      
    }
    nscore += jmax/(A.irow[i+1] -A.irow[i]); 
  }
  nscore /= double(n);
  std::cout<<" n score: "<<nscore <<std::endl;
  */

  Vector x(n, 0.0);
  Vector b(n, 1.0);
  /*  std::cout<<" start computing... " <<std::endl;
  for(int iter=0;  iter<100; ++iter){
    A.mult(x, b); 
  }
  */
  /*JacobiPreconditioner prec(A) ; 
  int n = A.n ;
  Vector x(n, 0.0);
  Vector b(n, 1.0);
  int nb_pcg_iter = pcg(A, prec, b, x, 5000, 1.e-16); 
  std::cout<<" nb iter: "<< nb_pcg_iter << std::endl ; 
  Vector res(n, 0.0);
  A.mult(x, res);
  axpy(res, -1, b); 
  std::cout<<" residual norm: "<< norm2(res) << std::endl;
  */
  /*// 1D Laplacian
  int n = 1000 ;
  Matrix A ;
  A.n = n ;
  A.nnz = 3*n-2;
  A.v = new double[3*n-2] ; 
  A.irow = new int[n+1] ;
  A.jcol = new int[3*n-2];
  int nnz = 0 ; 
  Vector b(n);
  for(int i=0; i<n; ++i){
    A.irow[i] = nnz ; 
    b[i] = 1.0;
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

  
  JacobiPreconditioner prec(A) ; 
  Vector x(n, 0.0);
  int nb_pcg_iter = pcg(A, prec, b, x, 1000, 1.e-6); 
  std::cout<<" nb iter: "<< nb_pcg_iter << std::endl ; 
  Vector res(n, 0.0);
  A.mult(x, res);
  std::cout<<" res: "<< res << std::endl;
  */
}; 
