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
    void set_value(const int idx, const fpctype value){gf[idx] = value;};
    fptype get_time(const int idx){return timepoints[idx];};
    fpctype get_value(const int idx){return gf[idx];};
    fpctype get_interpolated_value(const fptype tau);
    int get_ntimes(){return ntimes;};
  private:
    vector<fptype> timepoints;
    vector<fpctype> gf;
    int ntimes;
    fptype beta, dtau;
};

#endif