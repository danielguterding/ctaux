//files.hpp
#include <string>
#include <vector>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "typedefs.hpp"

using namespace std;

#ifndef __FILES_H_INCLUDED__
#define __FILES_H_INCLUDED__

struct SolverParameters{
  fptype beta;
  int nsamples;
};

class SolverParametersReader{
  public:
    SolverParametersReader(string infilename);
    SolverParameters get_parameters(){return p;};
  private:
    void read();
    string infilename;
    SolverParameters p;
};

#endif 

string trim_all(const std::string &str);