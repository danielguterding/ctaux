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
    ImaginaryTimeGreensFunction weissfield;
    gfreader.read_gf(p.inputfilepathweiss, weissfield);
    cout << "Input Weiss field successfully read." << endl;
    
    ImaginaryTimeGreensFunction outputgf;
    CTAUXSolver solver(p, weissfield, outputgf);
    cout << "Impurity solver successfully started." << endl;
    
    
    
  }
  else{
    cout << "Wrong number of input arguments.\nUsage: ctaux [string solverparameterinfilename]" << endl;
  }
  
  return 0;
}
