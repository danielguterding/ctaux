//random.hpp
#include <random>
#include "typedefs.hpp"

using namespace std;

// class RandomNumberGenerator{
//   public:
//     RandomNumberGenerator(const int key){}; //save an identifier that is used to initialize the RNG, e.g. node rank
//     virtual fptype get_value(){return 0;}; //this function returns a floating point value between 0 and 1
// };

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