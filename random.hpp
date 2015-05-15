//random.hpp
#include <random>

#include "typedefs.hpp"

using namespace std;

#ifndef __RANDOM_H_INCLUDED__
#define __RANDOM_H_INCLUDED__

class RNG_StdMersenne{
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

#endif