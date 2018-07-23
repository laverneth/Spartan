// -*- C++ -*-
#include <fstream>


namespace spartan_utils {

  static void load(const std::string& filename) {	
    std::ifstream myfile (filename);
    if (myfile.is_open()){
      int nrow, ncol, nnz ;
      myfile>>nrow>>ncol>>nnz ; 
      std::cout<<" nrow: "<< nrow << std::endl ;
      std::cout<<" nnz: "<< nnz << std::endl ;
      myfile.close();
    } 
  }
  

}
