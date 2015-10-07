//main.cpp
#include <boost/mpi.hpp>

#include "gf.hpp"
#include "solver.hpp"
#include "random.hpp"
#include "typedefs.hpp"
#include "files.hpp"

using namespace std;
namespace mpi = boost::mpi;

template <typename T>
struct elementwise_add {
    std::vector<T> operator()(const std::vector<T>& a, const std::vector<T>& b) const {
        std::vector<T> result(a.size());
        std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::plus<T>());
        return result;
    }
};

int main(int argc, char* argv[]){
  
  mpi::environment mpienv(argc, argv);
  mpi::communicator mpicomm;
  uint nnodes = mpicomm.size();
  uint noderank = mpicomm.rank();
  
  if(2 == argc){
    string solverparameterinfilename = argv[1];
    
    SolverParametersReader sreader(solverparameterinfilename);
    SolverParameters p = sreader.get_parameters();
    //distribute number of measurement steps across all nodes
    p.nsamplesmeasure = max(p.nsamplesmeasure/(long int) nnodes, (long int) 1);
    if(mpicomm.rank() == 0){cout << "Solver parameters successfully read." << endl;};
    
    ImaginaryTimeGreensFunctionReader gfreader;
    ImaginaryTimeGreensFunction weissfield_up, weissfield_dn;
    gfreader.read_gf(p.inputfilepathweiss_up, weissfield_up);
    gfreader.read_gf(p.inputfilepathweiss_dn, weissfield_dn);
    if(mpicomm.rank() == 0){cout << "Input Weiss field successfully read." << endl;};
    
    ImaginaryTimeGreensFunction outputgf_up, outputgf_dn;
    CTAUXSolver solver(p, weissfield_up, weissfield_dn, outputgf_up, outputgf_dn, noderank);
    if(mpicomm.rank() == 0){cout << "Impurity solver successfully started." << endl;};
    
    solver.do_warmup();
    if(mpicomm.rank() == 0){cout << "Warmup phase completed." << endl;};
    
    solver.do_measurement();
    if(mpicomm.rank() == 0){cout << "Measurement phase completed." << endl;};
    solver.construct_interacting_gf();
    
    //collect Greens function data
    vector<fptype> gflocal_up = outputgf_up.get_gf();
    vector<fptype> gflocal_dn = outputgf_dn.get_gf();
    std::transform(gflocal_up.begin(), gflocal_up.end(), gflocal_up.begin(), std::bind2nd(std::divides<fptype>(), fptype(nnodes))); //divide by number of nodes before gather operation to ensure averaging
    std::transform(gflocal_dn.begin(), gflocal_dn.end(), gflocal_dn.begin(), std::bind2nd(std::divides<fptype>(), fptype(nnodes)));    
    vector<fptype> gfaveraged_up(0), gfaveraged_dn(0);
    mpi::reduce(mpicomm, gflocal_up, gfaveraged_up, elementwise_add<fptype>(), 0);
    mpi::reduce(mpicomm, gflocal_dn, gfaveraged_dn, elementwise_add<fptype>(), 0);    
    mpicomm.barrier();
    if(mpicomm.rank() == 0){
      outputgf_up.set_gf(gfaveraged_up);
      outputgf_dn.set_gf(gfaveraged_dn);
    }
    
    if(mpicomm.rank() == 0){cout << "Constructed interacting Green's function from binned data." << endl;};
    if(mpicomm.rank() == 0){cout << "Impurity solver successfully finished." << endl;};
    
    //calculate averaged perturbation order additionally averaged over nodes
    const fptype polocal = solver.get_average_perturbation_order();
    fptype poaveraged = 0;
    mpi::reduce(mpicomm, polocal/nnodes, poaveraged, std::plus<fptype>(), 0);
    if(mpicomm.rank() == 0){cout << "Average perturbation order: " << poaveraged << endl;};
    
    //calculate averaged acceptance ratio over nodes
    const fptype arlocal = solver.get_acceptance_ratio();
    fptype araveraged = 0;
    mpi::reduce(mpicomm, arlocal/nnodes, araveraged, std::plus<fptype>(), 0);
    if(mpicomm.rank() == 0){cout << "Acceptance ratio: " << araveraged << endl;};
    
    if(mpicomm.rank() == 0){
      ImaginaryTimeGreensFunctionWriter gfwriter;
      gfwriter.write_gf(p.outputfilepathgf_up, outputgf_up);
      gfwriter.write_gf(p.outputfilepathgf_dn, outputgf_dn);
      cout << "Output Green's function successfully written." << endl;
      cout << "Program finished." << endl;
    }
    mpicomm.barrier();
  }
  else{
    cout << "Wrong number of input arguments.\nUsage: ctaux [string solverparameterinfilename]" << endl;
  }
  
  return 0;
}
