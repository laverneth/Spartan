#include "hsc_preconditioner.h"

namespace spartan {

  HSCPreconditioner::HSCPreconditioner(Matrix& A) {
    set_default_parameters(); 
    int level = 1 ;

     struct HierarchyLevel {
      int cycle ;
      unsigned int *C;
      Matrix P;
      Matrix A;
      double *invD;
      int fine_n, coarse_n;
      Matrix R;
      Matrix Rt;
      double *e, *r, *Ae;
      Cholesky chol;
    }; 


     /*while(level<= max_levels && hierarchy[level].size()>coarsest_level){


    }
    */
  }

  void HSCPreconditioner::apply(const Vector& x, Vector& y) const {

  }

  void HSCPreconditioner::set_default_parameters(){
    coarsest_level = 1024 ;
    pre_smooth_iter = 1;
    post_smooth_iter = 1;
    cycle = 1 ;
    max_levels = 1000; 
  }
  
  void HSCPreconditioner::set_param(const std::string& param_name, int value){
    if(param_name == "coarsest_level"){
      coarsest_level = value;
    }
    else if("pre_smooth_iter"){
      pre_smooth_iter = value;
    }
    else if("post_smooth_iter"){
      post_smooth_iter=value;
    }  
    else if("cycle"){
      cycle=value;
    }  
    else if("max_levels"){
      max_levels=value;
    }
  }

}
