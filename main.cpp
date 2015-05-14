//main.cpp
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
    
    
  }
  else{
    cout << "Wrong number of input arguments.\nUsage: ctaux [string solverparameterinfilename]" << endl;
  }
  
  return 0;
}
