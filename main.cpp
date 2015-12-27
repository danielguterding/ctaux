//main.cpp
#include <mpi.h>

#include "gf.hpp"
#include "solver.hpp"
#include "random.hpp"
#include "typedefs.hpp"
#include "files.hpp"

using namespace std;

int main(int argc, char* argv[]){
  
  MPI::Init(argc, argv);
  const uint nnodes = MPI::COMM_WORLD.Get_size();
  const uint noderank = MPI::COMM_WORLD.Get_rank();
  
  if(2 == argc){
    string solverparameterinfilename = argv[1];
    
    SolverParametersReader sreader(solverparameterinfilename);
    SolverParameters p = sreader.get_parameters();
    //distribute number of measurement steps across all nodes
    p.nsamplesmeasure = max(p.nsamplesmeasure/(long int) nnodes, (long int) 1);
    if(noderank == 0){cout << "Solver parameters successfully read." << endl;};
    
    ImaginaryTimeGreensFunctionReader gfreader;
    ImaginaryTimeGreensFunction weissfield_up, weissfield_dn;
    gfreader.read_gf(p.inputfilepathweiss_up, weissfield_up);
    gfreader.read_gf(p.inputfilepathweiss_dn, weissfield_dn);
    if(noderank == 0){cout << "Input Weiss field successfully read." << endl;};
    
    ImaginaryTimeGreensFunction outputgf_up, outputgf_dn;
    LegendreCoefficientRepresentation outputgf_legendre_up, outputgf_legendre_dn;
    CTAUXSolver solver(p, weissfield_up, weissfield_dn, outputgf_up, outputgf_dn, outputgf_legendre_up, outputgf_legendre_dn, noderank);
    if(noderank == 0){cout << "Impurity solver successfully started." << endl;};
    
    solver.do_warmup();
    if(noderank == 0){cout << "Warmup phase completed." << endl;};
    
    solver.do_measurement();
    if(noderank == 0){cout << "Measurement phase completed." << endl;};
    solver.construct_interacting_gf();
    
    //collect Greens function data
    vector<fptype> gflocal_up = outputgf_up.get_gf();
    vector<fptype> gflocal_dn = outputgf_dn.get_gf();
    vector<fptype> legendrelocal_up = outputgf_legendre_up.get_coefficients();
    vector<fptype> legendrelocal_dn = outputgf_legendre_dn.get_coefficients();
    std::transform(gflocal_up.begin(), gflocal_up.end(), gflocal_up.begin(), std::bind2nd(std::divides<fptype>(), fptype(nnodes))); //divide by number of nodes before gather operation to ensure averaging
    std::transform(gflocal_dn.begin(), gflocal_dn.end(), gflocal_dn.begin(), std::bind2nd(std::divides<fptype>(), fptype(nnodes))); 
    std::transform(legendrelocal_up.begin(), legendrelocal_up.end(), legendrelocal_up.begin(), std::bind2nd(std::divides<fptype>(), fptype(nnodes)));
    std::transform(legendrelocal_dn.begin(), legendrelocal_dn.end(), legendrelocal_dn.begin(), std::bind2nd(std::divides<fptype>(), fptype(nnodes)));
    vector<fptype> gfaveraged_up(gflocal_up.size(), 0), gfaveraged_dn(gflocal_dn.size(), 0), legendreaveraged_up(legendrelocal_up.size(), 0), legendreaveraged_dn(legendrelocal_dn.size(), 0);
    MPI::COMM_WORLD.Reduce(&gflocal_up.front(), &gfaveraged_up.front(), gflocal_up.size(), MPI_DOUBLE, MPI_SUM, 0);
    MPI::COMM_WORLD.Reduce(&gflocal_dn.front(), &gfaveraged_dn.front(), gflocal_dn.size(), MPI_DOUBLE, MPI_SUM, 0);
    MPI::COMM_WORLD.Reduce(&legendrelocal_up.front(), &legendreaveraged_up.front(), legendrelocal_up.size(), MPI_DOUBLE, MPI_SUM, 0);
    MPI::COMM_WORLD.Reduce(&legendrelocal_dn.front(), &legendreaveraged_dn.front(), legendrelocal_dn.size(), MPI_DOUBLE, MPI_SUM, 0);
    MPI::COMM_WORLD.Barrier();
    if(noderank == 0){
      outputgf_up.set_gf(gfaveraged_up);
      outputgf_dn.set_gf(gfaveraged_dn);
      outputgf_legendre_up.set_coefficients(legendreaveraged_up);
      outputgf_legendre_dn.set_coefficients(legendreaveraged_dn);
    }
    
    if(noderank == 0){cout << "Constructed interacting Green's function from binned data." << endl;};
    if(noderank == 0){cout << "Impurity solver successfully finished." << endl;};
    
    //calculate averaged perturbation order additionally averaged over nodes
    const fptype polocal = solver.get_average_perturbation_order();
    fptype poaveraged = 0;
    const fptype polocaltemp = polocal/nnodes;
    MPI::COMM_WORLD.Reduce(&polocaltemp, &poaveraged, 1, MPI_DOUBLE, MPI_SUM, 0);
    if(noderank == 0){cout << "Average perturbation order: " << poaveraged << endl;};
    
    //calculate averaged acceptance ratio over nodes
    const fptype arlocal = solver.get_acceptance_ratio();
    fptype araveraged = 0;
    const fptype arlocaltemp = arlocal/nnodes;
    MPI::COMM_WORLD.Reduce(&arlocaltemp, &araveraged, 1, MPI_DOUBLE, MPI_SUM, 0);
    if(noderank == 0){cout << "Acceptance ratio: " << araveraged << endl;};
    
    if(noderank == 0){
      ImaginaryTimeGreensFunctionWriter gfwriter;
      gfwriter.write_gf(p.outputfilepathgf_up, outputgf_up);
      gfwriter.write_gf(p.outputfilepathgf_dn, outputgf_dn);
      LegendreCoefficientRepresentationWriter coeffwriter;
      coeffwriter.write_coefficients(p.outputfilepathgf_up_legendre, outputgf_legendre_up);
      coeffwriter.write_coefficients(p.outputfilepathgf_dn_legendre, outputgf_legendre_dn);
      cout << "Output Green's function successfully written." << endl;
      cout << "Program finished." << endl;
    }
    MPI::COMM_WORLD.Barrier();
  }
  else{
    if(noderank == 0){
      cout << "Wrong number of input arguments.\nUsage: ctaux [string solverparameterinfilename]" << endl;
    }
  }
  
  MPI::Finalize();
  return 0;
}
