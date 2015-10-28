//gf.hpp
#include <vector>

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

#endif