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