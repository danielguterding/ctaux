//gf.cpp
#include "gf.hpp"

ImaginaryTimeGreensFunction::ImaginaryTimeGreensFunction(){
 
}

void ImaginaryTimeGreensFunction::generate_timeaxis(const fptype beta, const int ntimes){
  
  this->beta = beta;
  this->ntimes = ntimes;
  this->dtau = this->beta/(this->ntimes - 1);
  this->timepoints.resize(this->ntimes);
  this->gf.resize(this->ntimes);
  for(int i=0;i<this->ntimes;i++){
    timepoints[i] = i*dtau;
    gf[i] = 0;
  }
}

fptype ImaginaryTimeGreensFunction::get_interpolated_value(const fptype tau){
  //obtains values between the actually stored values of the gf by interpolating stored values linearly
  //antiperiodicity of the GF in imaginary time is automatically taken care of
  int p = floor(tau/beta); //counts the number of periods we have to translate back 
  fptype tauprime = tau - p*beta; //this is the timepoint in the original interval where we want to look up value 
  int lidx = floor(tauprime/dtau); //this is the lower index of the stored timepoints neighbouring the timepoint of interest
  fptype taurem = tauprime - lidx*dtau; //this is the remainder of tauprime within the interval of interest
  fptype gfval = gf[lidx] + (gf[lidx+1] - gf[lidx])/dtau*taurem;
  fptype sign = pow(-1,p); //get the sign corresponding to fermionic antiperiodicity
  return sign*gfval;
}