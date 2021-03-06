//gf.hpp
#include <vector>
#include <boost/math/special_functions/legendre.hpp>

#include "typedefs.hpp"

using namespace std;

#ifndef __GREENSFUNCTION_H_INCLUDED__
#define __GREENSFUNCTION_H_INCLUDED__

class ImaginaryTimeGreensFunction{
  public:
    ImaginaryTimeGreensFunction();
    void generate_timeaxis(const fptype beta, const int ntimes);
    void set_value(const int idx, const fptype value){gf[idx] = value;};
    fptype get_time(const int idx){return timepoints[idx];};
    fptype get_value(const int idx){return gf[idx];};
    fptype get_interpolated_value(const fptype tau);
    int get_ntimes(){return ntimes;};
    vector<fptype> get_gf(){return gf;};
    void set_gf(vector<fptype> gf){this->gf = gf;};
  private:
    vector<fptype> timepoints;
    vector<fptype> gf;
    int ntimes;
    fptype beta, dtau;
};

class LegendreCoefficientRepresentation{
  public:
    LegendreCoefficientRepresentation();
    void initialize(const int ncoeff, const fptype beta);
    fptype get_coefficient(const int idx){return coefficients[idx];};
    void set_coefficient(const int idx, const fptype value){coefficients[idx] = value;};
    fptype legendre_p(const int order, const fptype x);
    vector<fptype> get_coefficients(){return coefficients;};
    void set_coefficients(const vector<fptype> coefficients){this->coefficients = coefficients;};
    int get_ncoefficients(){return ncoeff;};
    fptype x(const fptype tau){return 2*tau/beta - 1;};
    fptype t(const int l, const int p);
    void calculate_imaginary_time_gf_from_sigmag(ImaginaryTimeGreensFunction& g0, ImaginaryTimeGreensFunction& sigmag);
  private:
    int ncoeff;
    fptype beta;
    vector<fptype> coefficients;
    long long int factorial(long long int n){return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;};
};

#endif