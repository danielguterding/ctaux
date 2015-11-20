//files.hpp
#include <string>
#include <vector>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "gf.hpp"
#include "typedefs.hpp"

using namespace std;

#ifndef __FILES_H_INCLUDED__
#define __FILES_H_INCLUDED__

struct SolverParameters{
  fptype beta, K, U;
  long int nsampleswarmup, nsamplesmeasure, nbins;
  int nlegendre;
  string inputfilepathweiss_up, inputfilepathweiss_dn, outputfilepathgf_up, outputfilepathgf_dn, outputfilepathgf_up_legendre, outputfilepathgf_dn_legendre;
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

class ImaginaryTimeGreensFunctionReader{
  public:
    ImaginaryTimeGreensFunctionReader();
    void read_gf(const string infilename, ImaginaryTimeGreensFunction& gf);
  private:
    
};

class ImaginaryTimeGreensFunctionWriter{
  public:
    ImaginaryTimeGreensFunctionWriter();
    void write_gf(const string outfilename, ImaginaryTimeGreensFunction& gf);
  private:
    
};

class LegendreCoefficientRepresentationWriter{
  public:
    LegendreCoefficientRepresentationWriter();
    void write_coefficients(string outfilename, LegendreCoefficientRepresentation& legendre);
  private:
  
};

#endif 

string trim_all(const std::string &str);