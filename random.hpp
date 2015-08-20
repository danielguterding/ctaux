//random.hpp
#include <random>
#include <Random123/philox.h>

#include "typedefs.hpp"

using namespace std;

#ifndef __RANDOM_H_INCLUDED__
#define __RANDOM_H_INCLUDED__

class RNG{ //base class that contains the virtual public method that all derived classes have to implement
  public:
    virtual ~RNG(){};
    virtual fptype get_value(){return 0;};
};

class RNG_StdMersenne : public RNG{
  //Mersenne Twister from the C++ standard library
  public:
    RNG_StdMersenne(const int key);
    ~RNG_StdMersenne();
    fptype get_value();
  private:
    int key;
    std::mt19937 *rng_ptr;
    std::uniform_real_distribution<fptype> *dist_ptr;
};

class RNG_Philox4x32 : public RNG{
  //Philox4x32 from Random123 library
  public:
    RNG_Philox4x32(const int key);
    ~RNG_Philox4x32();
    fptype get_value();
  private:
    int valcounter;
    long int totalcounter;
    fptype maxval;
    fptype returnval;
    r123::Philox4x32* rng_ptr;
    r123::Philox4x32::ctr_type c={{}}; //contains loop variable
    r123::Philox4x32::ctr_type r={{}}; //contains random numbers
    r123::Philox4x32::key_type k={{}}; //contains seed
};

#endif