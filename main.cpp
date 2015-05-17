//main.cpp
#include "gf.hpp"
#include "solver.hpp"
#include "random.hpp"
#include "typedefs.hpp"
#include "files.hpp"

using namespace std;

int main(int argc, char* argv[]){
  
  if(2 == argc){
    string solverparameterinfilename = argv[1];
    
    SolverParametersReader sreader(solverparameterinfilename);
    SolverParameters p = sreader.get_parameters();
    cout << "Solver parameters successfully read." << endl;
    
    ImaginaryTimeGreensFunctionReader gfreader;
    ImaginaryTimeGreensFunction weissfield_up, weissfield_dn;
    gfreader.read_gf(p.inputfilepathweiss_up, weissfield_up);
    gfreader.read_gf(p.inputfilepathweiss_dn, weissfield_dn);
    cout << "Input Weiss field successfully read." << endl;
    
    ImaginaryTimeGreensFunction outputgf_up, outputgf_dn;
    CTAUXSolver solver(p, weissfield_up, weissfield_dn, outputgf_up, outputgf_dn);
    cout << "Impurity solver successfully started." << endl;
    
    solver.do_warmup();
    cout << "Warmup phase completed." << endl;
    
    solver.do_measurement();
    cout << "Measurement phase completed." << endl;
    solver.construct_interacting_gf();
    cout << "Constructed interacting Green's function from binned data." << endl;
    cout << "Impurity solver successfully finished." << endl;
    cout << "Average perturbation order: " << solver.get_average_perturbation_order() << endl;
    
    ImaginaryTimeGreensFunctionWriter gfwriter;
    gfwriter.write_gf(p.outputfilepathgf_up, outputgf_up);
    gfwriter.write_gf(p.outputfilepathgf_dn, outputgf_dn);
    cout << "Output Green's function successfully written." << endl;
    cout << "Program finished." << endl;
  }
  else{
    cout << "Wrong number of input arguments.\nUsage: ctaux [string solverparameterinfilename]" << endl;
  }
  
  return 0;
}
