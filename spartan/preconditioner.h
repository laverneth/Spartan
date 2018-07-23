// -*- C++ -*-
//

#ifndef __SPARTAN_preconditioner_h__
#define __SPARTAN_preconditioner_h__

#include "vector.h"
#include "sparse_types.h"

namespace spartan {

  /**
   * @brief Base class for preconditioners
   */
  class Preconditioner {
  public:
    ~Preconditioner(){};
    virtual void apply(const Vector& x, Vector& y) const = 0 ; 
  }; 
  

  class DiagonalPreconditioner: public Preconditioner {
  public:
    DiagonalPreconditioner(const Matrix& A); 
    ~DiagonalPreconditioner(){}
    virtual void apply(const Vector& x, Vector& y) const ; 
  private:
    Vector diag_ ; 
  }; 


}
#endif
