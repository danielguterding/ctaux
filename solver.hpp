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
    CTAUXSolver(SolverParameters& p, ImaginaryTimeGreensFunction& weissfield_up, ImaginaryTimeGreensFunction& weissfield_dn, ImaginaryTimeGreensFunction& outputgf_up, ImaginaryTimeGreensFunction& outputgf_dn, LegendreCoefficientRepresentation& outputgf_legendre_up, LegendreCoefficientRepresentation& outputgf_legendre_dn, const uint noderank);
    ~CTAUXSolver();
    void do_warmup();
    void do_measurement();
    void construct_interacting_gf();
    fptype get_average_perturbation_order(){return average_po;};
    fptype get_acceptance_ratio(){return acceptance_ratio;};
  private:
    void initialize();
    void step();
    void insert_update();
    void remove_update();
    void measure_gf();
    void measure_perturbation_order();
    void measure_acceptance_ratio();
    #if DEBUG 
    void calculate_Ninverse();
    #endif
    fptype egamma(int physicalspin, int auxiliaryspin);
    fptype gammaparameter;
    uint noderank;
    SolverParameters p;
    ImaginaryTimeGreensFunction *wfup_ptr, *wfdn_ptr, *outputgfup_ptr, *outputgfdn_ptr;
    LegendreCoefficientRepresentation *outputgflegendreup_ptr, *outputgflegendredn_ptr;
    CTAUXConfiguration *config_ptr;
    RNG *rng_ptr;
    Eigen::MatrixXd Nmatup, Nmatdn;
    fptype binwidth;
    vector<fptype> binmids;
    vector<fptype> gfupbins, gfdnbins, gfuplegendre, gfdnlegendre;
    fptype average_po, acceptance_ratio;
    bool update_accepted;
};

#endif