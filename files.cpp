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
  vector<string> splitline;
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.beta = boost::lexical_cast<fptype>(line);
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.nsamples = int(boost::lexical_cast<fptype>(line));
}

string trim_all(const std::string &str){  //with a more recent version of boost boost::trim_all() can be used instead of this function

  return boost::algorithm::find_format_all_copy(
    boost::trim_copy(str),
    boost::algorithm::token_finder (boost::is_space(),boost::algorithm::token_compress_on),
    boost::algorithm::const_formatter(" "));
}
