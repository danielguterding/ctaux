//solver.hpp
#include "gf.hpp"
#include "random.hpp"
#include "files.hpp"

#ifndef __SOLVER_H_INCLUDED__
#define __SOLVER_H_INCLUDED__

class CTAUXSolver{
  public:
    CTAUXSolver(SolverParameters& p, ImaginaryTimeGreensFunction& weissfield, ImaginaryTimeGreensFunction& outputgf);
    ~CTAUXSolver();
  private:
    SolverParameters p;
    ImaginaryTimeGreensFunction *weissfield_ptr, *outputgf_ptr;
};

#endif