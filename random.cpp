//random.cpp
#include "random.hpp"

RNG_StdMersenne::RNG_StdMersenne(const int key){
  
  this->key = key;
  this->rng_ptr = new std::mt19937(this->key);
  this->dist_ptr = new std::uniform_real_distribution<fptype>(0,1);
}

RNG_StdMersenne::~RNG_StdMersenne(){
  
  delete rng_ptr;
  delete dist_ptr;
}

fptype RNG_StdMersenne::get_value(){
  
  return (*dist_ptr)(*rng_ptr);
}

RNG_Philox4x32::RNG_Philox4x32(const int key){
  
  this->rng_ptr = new r123::Philox4x32;
  this->k[0] = key;
  valcounter = 3; //Philox4x32 generates four 32bit random values in each round, this counter determines how many of the four have been used
  totalcounter = 0; //this variable is a running index of how often the random number generator has been executed
  maxval = pow(2,32)-1;
}

RNG_Philox4x32::~RNG_Philox4x32(){
  
  delete rng_ptr;
}

fptype RNG_Philox4x32::get_value(){
  if(4 == valcounter){ //regenerate random values if all four have been used
    c[0] = totalcounter;
    r = (*rng_ptr)(c, k);
    valcounter = 0;
    totalcounter++;
  }
  returnval = r[valcounter++]/maxval;
  return returnval;
}