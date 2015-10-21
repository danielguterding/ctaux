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
  if(fabs(tau - beta) < 1e-13){ //fix last point which is tau=beta to be picked from first perio regardless
    p = 0;
  }
  const fptype tauprime = tau - p*beta; //this is the timepoint in the original interval where we want to look up value 
  const int lidx = floor(tauprime/dtau); //this is the lower index of the stored timepoints neighbouring the timepoint of interest
  const fptype taurem = tauprime - lidx*dtau; //this is the remainder of tauprime within the interval of interest
  const fptype gfval = gf[lidx] + (gf[lidx+1] - gf[lidx])/dtau*taurem;
  const fptype sign = pow(-1,p); //get the sign corresponding to fermionic antiperiodicity
  return sign*gfval;
}

LegendreCoefficientRepresentation::LegendreCoefficientRepresentation(){
  
}

void LegendreCoefficientRepresentation::initialize(const int ncoeff, const fptype beta){
  
  this->ncoeff = ncoeff;
  this->beta = beta;
  this->coefficients = vector<fptype>(this->ncoeff, 0);
}

fptype LegendreCoefficientRepresentation::legendre_p(const int order, const fptype x){
  
  return boost::math::legendre_p<fptype>(order, x);
}

fptype LegendreCoefficientRepresentation::t(const int l, const int p){
  
  if(0 == (p+l)%2){ //for p+l even return zero
    return 0;
  }
  else if(1 == p){ //use analytic simplification for case 1==p
    return -2*sqrt(2.0*l+1.0);
  }
  else{
    return pow(-1,p)*2*sqrt(2*l+1)*factorial(l+p-1)/(factorial(p-1)*factorial(l-p+1));
  }
}

void LegendreCoefficientRepresentation::calculate_imaginary_time_gf_from_sigmag(ImaginaryTimeGreensFunction& g0, ImaginaryTimeGreensFunction& sigmag){
  
  const int ntimes = sigmag.get_ntimes();
  for(int i=0;i<ntimes;i++){
    const fptype tau = sigmag.get_time(i);
    const fptype xval = this->x(tau);
    fptype gfval = 0;
    const fptype g0val = g0.get_interpolated_value(tau);
    //cout << g0val << endl;
    for(int j=0;j<ncoeff;j++){
      gfval += sqrt(2.0*j+1.0)/beta*legendre_p(j, xval)*coefficients[j];
    }
    gfval = g0val*(1.0+gfval);
    sigmag.set_value(i, gfval);
  }
}