//solver.hpp
#include "random.hpp"
#include "files.hpp"

#ifndef __SOLVER_H_INCLUDED__
#define __SOLVER_H_INCLUDED__

class CTAUXSolver{
  public:
    CTAUXSolver(SolverParameters& p);
    ~CTAUXSolver();
  private:
    SolverParameters p;
};

#endif