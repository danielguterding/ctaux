//files.cpp
#include "files.hpp"

SolverParametersReader::SolverParametersReader(string infilename){
  
  this->infilename = infilename;
  read();
}

void SolverParametersReader::read(){
  
  boost::filesystem::path infilepath(infilename);
  if(!boost::filesystem::exists(infilepath)){cout << "Warning! Input file does not exist!" << endl;};
  boost::filesystem::ifstream infilehandle(infilepath);
  
  string line = "";
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.beta = boost::lexical_cast<fptype>(line);
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.K = boost::lexical_cast<fptype>(line);
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.U = boost::lexical_cast<fptype>(line);
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.nsampleswarmup = int(boost::lexical_cast<fptype>(line));
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.nsamplesmeasure = int(boost::lexical_cast<fptype>(line));
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.inputfilepathweiss = line;
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.outputfilepathgf = line;
  
  infilehandle.close();
}

ImaginaryTimeGreensFunctionReader::ImaginaryTimeGreensFunctionReader(){
  
}

void ImaginaryTimeGreensFunctionReader::read_gf(const string infilename, ImaginaryTimeGreensFunction& gf){
  //first load the input file
  boost::filesystem::path infilepath(infilename);
  if(!boost::filesystem::exists(infilepath)){cout << "Warning! Input file does not exist!" << endl;};
  boost::filesystem::ifstream infilehandle(infilepath);
  
  string line = "";
  vector<string> splitline;
  
  vector<fptype> taupoints;
  vector<fpctype> inputgf;
  
  getline(infilehandle, line); //read comment line and discard it
  while(getline(infilehandle, line)){
    line = trim_all(line); //erase trailing and leading spaces, reduce intermediate spaces to one space
    boost::split(splitline, line, boost::is_any_of("\t "));
    taupoints.push_back(boost::lexical_cast<fptype>(splitline[0]));
    inputgf.push_back(fpctype(boost::lexical_cast<fptype>(splitline[1]), boost::lexical_cast<fptype>(splitline[2])));
  }
  infilehandle.close();
  
  const int ntimes = taupoints.size();
  const fptype beta = taupoints[ntimes-1];
  
  gf.generate_timeaxis(beta, ntimes);
  for(int i=0;i<ntimes;i++){
    gf.set_value(i, inputgf[i]);
  }
}

string trim_all(const std::string &str){  //with a more recent version of boost boost::trim_all() can be used instead of this function

  return boost::algorithm::find_format_all_copy(
    boost::trim_copy(str),
    boost::algorithm::token_finder (boost::is_space(),boost::algorithm::token_compress_on),
    boost::algorithm::const_formatter(" "));
}
