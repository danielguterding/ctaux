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
  this->p.nsampleswarmup = (long int)(boost::lexical_cast<fptype>(line));
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.nsamplesmeasure = (long int)(boost::lexical_cast<fptype>(line));
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.nbins = (long int)(boost::lexical_cast<fptype>(line));
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.nlegendre= (int)(boost::lexical_cast<fptype>(line));
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.inputfilepathweiss_up = line;
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.inputfilepathweiss_dn = line;
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.outputfilepathgf_up = line;
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.outputfilepathgf_dn = line;
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.outputfilepathgf_up_legendre = line;
  
  getline(infilehandle, line); //read comment line and discard it
  getline(infilehandle, line);
  this->p.outputfilepathgf_dn_legendre = line;
  
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
  vector<fptype> inputgf;
  
  getline(infilehandle, line); //read comment line and discard it
  while(getline(infilehandle, line)){
    line = trim_all(line); //erase trailing and leading spaces, reduce intermediate spaces to one space
    boost::split(splitline, line, boost::is_any_of("\t "));
    taupoints.push_back(boost::lexical_cast<fptype>(splitline[0]));
    //input GF on tau axis must be positive!
    inputgf.push_back(fabs(boost::lexical_cast<fptype>(splitline[1])));
  }
  infilehandle.close();
  
  const int ntimes = taupoints.size();
  const fptype beta = taupoints[ntimes-1];
  
  gf.generate_timeaxis(beta, ntimes);
  for(int i=0;i<ntimes;i++){
    gf.set_value(i, inputgf[i]);
  }
}

ImaginaryTimeGreensFunctionWriter::ImaginaryTimeGreensFunctionWriter(){
  
}

void ImaginaryTimeGreensFunctionWriter::write_gf(const string outfilename, ImaginaryTimeGreensFunction& gf){
  
  boost::filesystem::path outfilepath(outfilename);
  boost::filesystem::ofstream outfilehandle(outfilepath);
  
  const int ntimes = gf.get_ntimes();
  
  outfilehandle << "#tau in eV^-1 from 0 to beta; GF in eV^-1 real part, imag part is zero;\n";
  fptype val = 0;
  for(int i=0;i<ntimes;i++){
    val = -fabs(gf.get_value(i));
    outfilehandle << boost::lexical_cast<string>(boost::format("%1.14e % 1.14e\n") % gf.get_time(i) % val);
  }
  
  outfilehandle.close();
}

LegendreCoefficientRepresentationWriter::LegendreCoefficientRepresentationWriter(){
  
}

void LegendreCoefficientRepresentationWriter::write_coefficients(string outfilename, LegendreCoefficientRepresentation& legendre){
  
  boost::filesystem::path outfilepath(outfilename);
  boost::filesystem::ofstream outfilehandle(outfilepath);
  
  outfilehandle << "#idx, coefficient value\n";
  for(int i=0;i<legendre.get_ncoefficients();i++){
    outfilehandle << boost::lexical_cast<string>(boost::format("%03i % 1.14e\n") % i % legendre.get_coefficient(i));
  }
  
  outfilehandle.close();
}

string trim_all(const std::string &str){  //with a more recent version of boost boost::trim_all() can be used instead of this function

  return boost::algorithm::find_format_all_copy(
    boost::trim_copy(str),
    boost::algorithm::token_finder (boost::is_space(),boost::algorithm::token_compress_on),
    boost::algorithm::const_formatter(" "));
}
