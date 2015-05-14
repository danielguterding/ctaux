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
    
    
//     RNG_StdMersenne rng(10);
//     for(int i=0;i<100;i++){
//       cout << rng.get_value() << endl;
//     }
  }
  else{
    cout << "Wrong number of input arguments.\nUsage: ctaux [string solverparameterinfilename]" << endl;
  }
  
  return 0;
}
