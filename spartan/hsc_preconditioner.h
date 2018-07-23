// -*- C++ -*-
//

#ifndef __SPARTAN_hsc_preconditioner_h__
#define __SPARTAN_hsc_preconditioner_h__

#include "preconditioner.h"
#include "sparse_types.h"

namespace spartan {

  /**
   * @brief Hierarchical Sparsify Compensate preconditioner
   */
  class HSCPreconditioner : public Preconditioner {
  public:
    HSCPreconditioner(Matrix& A) ; 
    virtual void apply(const Vector& x, Vector& y) const ;

    void set_default_parameters(); 
    void set_param(const std::string& param_name, int param_value); 

  private:
    int coarsest_level ;
    int pre_smooth_iter;
    int post_smooth_iter;
    int cycle;
    int max_levels; 

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

    HierarchyLevel* levels_ ; 
  }; 



}
#endif
