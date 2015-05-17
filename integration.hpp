//integration.hpp
#include <iostream>
#include <vector>

#include "typedefs.hpp"

using namespace std;

#ifndef __INTEGRATION_H_INCLUDED__
#define __INTEGRATION_H_INCLUDED__

template <class T> 
class TrapezoidalIntegrator{
  public:
    TrapezoidalIntegrator(){};
    void set_x(vector<fptype> x){this->x = x;};
    void set_y(vector<T> y){this->y = y;};
    void calculate(){
      
      value = 0;
      
      if(x.size() != y.size()){
	cout << "Size of integration array x does not match size of y array. Aborting." << endl;
      }
      else{
	for(uint i=0;i<x.size()-1;i++){
	  value += 0.5*(x[i+1] - x[i])*(y[i+1] + y[i]);
	}
      }
    };
    T get_value(){return value;};
  private:
    vector<fptype> x;
    vector<T> y;
    T value;
};

#endif