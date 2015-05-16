//solver.cpp
#include "solver.hpp"

CTAUXSolver::CTAUXSolver(SolverParameters& p, ImaginaryTimeGreensFunction& weissfield, ImaginaryTimeGreensFunction& outputgf){
  
  this->p = p;
  this->weissfield_ptr = &weissfield;
  this->outputgf_ptr = &outputgf;
  
  //initialize outputgf with beta and number of time slices, i.e. bins
}

CTAUXSolver::~CTAUXSolver(){
  
}