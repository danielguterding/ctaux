//solver.hpp
#include <Eigen/Dense>

#include "gf.hpp"
#include "random.hpp"
#include "files.hpp"

#ifndef __SOLVER_H_INCLUDED__
#define __SOLVER_H_INCLUDED__

class CTAUXConfiguration{
  public:
    CTAUXConfiguration();
    void insert(fptype tau, bool spin);
    void remove(const uint idx);
    int get_auxspin(const uint idx);
    fptype get_time(const uint idx);
    int get_perturbation_order(){return perturbationorder;};
  private:
    int perturbationorder;
    vector<bool> auxiliaryspins;
    vector<fptype> timepoints;
};

class CTAUXSolver{
  public:
    CTAUXSolver(SolverParameters& p, ImaginaryTimeGreensFunction& weissfield_up, ImaginaryTimeGreensFunction& weissfield_dn, ImaginaryTimeGreensFunction& outputgf_up, ImaginaryTimeGreensFunction& outputgf_dn);
    ~CTAUXSolver();
    void do_warmup();
    void do_measurement();
  private:
    void initialize();
    void step();
    void insert_update();
    void remove_update();
    fptype egamma(int physicalspin, int auxiliaryspin);
    fptype gammaparameter;
    SolverParameters p;
    ImaginaryTimeGreensFunction *wfup_ptr, *wfdn_ptr, *outputgfup_ptr, *outputgfdn_ptr;
    CTAUXConfiguration *config_ptr;
    RNG_StdMersenne *rng_ptr;
    Eigen::MatrixXcd Nmatup, Nmatdn;
};

#endif