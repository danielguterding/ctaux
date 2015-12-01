//solver.cpp
#include "solver.hpp"

CTAUXConfiguration::CTAUXConfiguration(){
  
  auxiliaryspins.resize(0);
  timepoints.resize(0);
  perturbationorder = 0;
}

void CTAUXConfiguration::insert(fptype tau, bool spin){
  //spin==0 means down, spin==1 means up
  timepoints.push_back(tau);
  auxiliaryspins.push_back(spin);
  perturbationorder++;
}

void CTAUXConfiguration::remove(const uint idx){
  
  timepoints.erase(timepoints.begin()+idx);
  auxiliaryspins.erase(auxiliaryspins.begin()+idx);
  perturbationorder--;
}

int CTAUXConfiguration::get_auxspin(const uint idx){
  //returns -1 for down spin, +1 for up spin
  if(auxiliaryspins[idx]){
    return 1;
  }
  else{
    return -1;
  }
}

fptype CTAUXConfiguration::get_time(const uint idx){
  
  return timepoints[idx];
}

CTAUXSolver::CTAUXSolver(SolverParameters& p, ImaginaryTimeGreensFunction& weissfield_up, ImaginaryTimeGreensFunction& weissfield_dn, ImaginaryTimeGreensFunction& outputgf_up, ImaginaryTimeGreensFunction& outputgf_dn, LegendreCoefficientRepresentation& outputgf_legendre_up, LegendreCoefficientRepresentation& outputgf_legendre_dn, const uint noderank){
  
  this->p = p;
  this->wfup_ptr = &weissfield_up;
  this->wfdn_ptr = &weissfield_dn;
  this->outputgfup_ptr = &outputgf_up;
  this->outputgfdn_ptr = &outputgf_dn;
  this->outputgflegendreup_ptr = &outputgf_legendre_up;
  this->outputgflegendredn_ptr = &outputgf_legendre_dn;
  
  this->noderank = noderank;
  
  //initialize outputgf with beta and number of time slices, i.e. bins
  this->outputgfup_ptr->generate_timeaxis(this->p.beta, this->p.nbins);
  this->outputgfdn_ptr->generate_timeaxis(this->p.beta, this->p.nbins);
  
  //initialize bins for conventional gf measurement
  this->binwidth = this->p.beta/fptype(this->p.nbins);
  this->binmids.resize(this->p.nbins);
  this->gfupbins.resize(this->p.nbins);
  this->gfdnbins.resize(this->p.nbins);
  for(int i=0;i<this->p.nbins;i++){
    this->binmids[i] = (0.5+i)*this->binwidth;
    this->gfupbins[i] = 0;
    this->gfdnbins[i] = 0;
  }
  
  //initialize Legendre coefficient bins
  this->gfuplegendre.resize(this->p.nlegendre);
  this->gfdnlegendre.resize(this->p.nlegendre);
  for(int i=0;i<this->p.nlegendre;i++){
    gfuplegendre[i] = 0;
    gfdnlegendre[i] = 0;
  }
  this->outputgflegendreup_ptr->initialize(this->p.nlegendre, this->p.beta);
  this->outputgflegendredn_ptr->initialize(this->p.nlegendre, this->p.beta);
  
  this->gammaparameter = acosh(1 + p.beta*p.U/2.0/p.K);
  this->update_accepted = 0;
  this->average_po = 0;
  this->acceptance_ratio = 0;
  
  initialize();
}

CTAUXSolver::~CTAUXSolver(){
  
  delete config_ptr;
  delete rng_ptr;
}

void CTAUXSolver::initialize(){
  
  const int seed = noderank; //use noderank as seed
  config_ptr = new CTAUXConfiguration;
  
  //choose random number generator
  const int rngidx = 1;
  switch(rngidx){
    case 0:
      rng_ptr = new RNG_StdMersenne(seed);
      break;
    case 1:
      rng_ptr = new RNG_Threefry4x64(seed);
      break;
    case 2:
      rng_ptr = new RNG_Philox4x64(seed);
      break;
    case 3:
      rng_ptr = new RNG_Philox4x32(seed);
      break;
  }
      
  const fptype tau = p.beta*rng_ptr->get_value();
  const bool spin = floor(2*rng_ptr->get_value());
  config_ptr->insert(tau, spin);
  const int auxspin = config_ptr->get_auxspin(0);
  
  //update N matrix for up and dn
  Nmatup.resize(1,1);
  Nmatdn.resize(1,1);
  Nmatup(0,0) = 1.0/(egamma( 1, auxspin) - (egamma( 1, auxspin) - 1)*wfup_ptr->get_interpolated_value(0));
  Nmatdn(0,0) = 1.0/(egamma(-1, auxspin) - (egamma(-1, auxspin) - 1)*wfdn_ptr->get_interpolated_value(0));
}

fptype CTAUXSolver::egamma(int physicalspin, int auxiliaryspin){
  
  return exp(gammaparameter*fptype(physicalspin)*fptype(auxiliaryspin));
}

void CTAUXSolver::do_warmup(){
  
  for(long int i=0;i<p.nsampleswarmup;i++){
    step();
  }
}

void CTAUXSolver::do_measurement(){
  
  for(long int i=0;i<p.nsamplesmeasure;i++){
    step();
    measure_gf();
    measure_perturbation_order();
    measure_acceptance_ratio();
    #if DEBUG
    calculate_Ninverse();
    #endif
  }
}

void CTAUXSolver::step(){
  
  update_accepted = 0;
  if(rng_ptr->get_value() < 0.5){
    insert_update();
  }
  else{
    remove_update();
  }
}

void CTAUXSolver::insert_update(){
  
  const fptype tau = p.beta*rng_ptr->get_value();
  const bool auxspin_dummy = floor(2*rng_ptr->get_value());
  const int auxspin = (auxspin_dummy ? 1 : -1); //if bool variable is 1, assign +1 spin, otherwise -1 spin
  const int po = config_ptr->get_perturbation_order();
  
  //calculate non-tilde quantities
  Eigen::VectorXd Qup = Eigen::VectorXd::Zero(po);
  Eigen::VectorXd Qdn = Eigen::VectorXd::Zero(po);
  Eigen::VectorXd Rup = Eigen::VectorXd::Zero(po);
  Eigen::VectorXd Rdn = Eigen::VectorXd::Zero(po);
  const fptype eup = egamma(1, auxspin);
  const fptype edn = egamma(-1, auxspin);
  const fptype Sup = eup - (eup - 1.0)*wfup_ptr->get_interpolated_value(0);
  const fptype Sdn = edn - (edn - 1.0)*wfdn_ptr->get_interpolated_value(0);
  for(int l=0;l<po;l++){
    const fptype taul = config_ptr->get_time(l);
    const int auxspinl = config_ptr->get_auxspin(l);
    Qup(l) = -(egamma( 1, auxspinl) - 1)*wfup_ptr->get_interpolated_value(taul - tau);
    Qdn(l) = -(egamma(-1, auxspinl) - 1)*wfdn_ptr->get_interpolated_value(taul - tau);
    Rup(l) = -(eup - 1.0)*wfup_ptr->get_interpolated_value(tau - taul);
    Rdn(l) = -(edn - 1.0)*wfdn_ptr->get_interpolated_value(tau - taul);
  }
  const fptype paccup = Sup - Rup.dot(Nmatup*Qup);
  const fptype paccdn = Sdn - Rdn.dot(Nmatdn*Qdn);
  const fptype pacc = min(1.0, p.K/(po+1.0)*fabs(paccup*paccdn));
  
  #if DEBUG
  const fptype detninvup_n = Nmatup.inverse().determinant();
  const fptype detninvdn_n = Nmatdn.inverse().determinant();
  #endif
  
  const fptype r = rng_ptr->get_value();
  if(r<pacc){ //accept update
    update_accepted = 1;
    config_ptr->insert(tau, auxspin_dummy);
    
    //calculate tilde quantities
    const fptype Stildeup = 1.0/paccup; 
    const fptype Stildedn = 1.0/paccdn;
    Eigen::VectorXd Qtildeup = -Nmatup*Qup*Stildeup;
    Eigen::VectorXd Qtildedn = -Nmatdn*Qdn*Stildedn;
    Eigen::VectorXd Rtildeup = -Stildeup*(Rup.transpose()*Nmatup).transpose();
    Eigen::VectorXd Rtildedn = -Stildedn*(Rdn.transpose()*Nmatdn).transpose();
    Eigen::MatrixXd Ptildeup = Nmatup + Qtildeup*Rtildeup.transpose()/Stildeup;
    Eigen::MatrixXd Ptildedn = Nmatdn + Qtildedn*Rtildedn.transpose()/Stildedn;
    Nmatup.noalias() = Eigen::MatrixXd::Zero(po+1, po+1);
    Nmatdn.noalias() = Eigen::MatrixXd::Zero(po+1, po+1);
    Nmatup.block(0,0,po,po).noalias() = Ptildeup;
    Nmatdn.block(0,0,po,po).noalias() = Ptildedn;
    Nmatup.block(0,po,po,1).noalias() = Qtildeup;
    Nmatdn.block(0,po,po,1).noalias() = Qtildedn;
    Nmatup.block(po,0,1,po).noalias() = Rtildeup.transpose();
    Nmatdn.block(po,0,1,po).noalias() = Rtildedn.transpose();
    Nmatup(po,po) = Stildeup;
    Nmatdn(po,po) = Stildedn;
    
    #if DEBUG
    const fptype detninvup_np1 = Nmatup.inverse().determinant();
    const fptype detninvdn_np1 = Nmatdn.inverse().determinant();
    const fptype pacc_explicit = min(1.0, p.K/(po+1.0)*detninvup_np1*detninvdn_np1/detninvup_n/detninvdn_n);
    cout << pacc << " " << pacc_explicit << endl;
    #endif
  }
}

void CTAUXSolver::remove_update(){
  
  const int po = config_ptr->get_perturbation_order();
  if(po > 0){
    const int ridx = floor(po*rng_ptr->get_value());
    //obtain tilde quantities directly from N matrices
    const fptype Stildeup = Nmatup(ridx, ridx);
    const fptype Stildedn = Nmatdn(ridx, ridx);
    Eigen::VectorXd Qtildeup = Eigen::VectorXd::Zero(po-1);
    Eigen::VectorXd Qtildedn = Eigen::VectorXd::Zero(po-1);
    Eigen::VectorXd Rtildeup = Eigen::VectorXd::Zero(po-1);
    Eigen::VectorXd Rtildedn = Eigen::VectorXd::Zero(po-1);
    for(int i=0;i<po;i++){
      if(i<ridx){
        Qtildeup(i) = Nmatup(i,ridx);
        Qtildedn(i) = Nmatdn(i,ridx);
        Rtildeup(i) = Nmatup(ridx,i);
        Rtildedn(i) = Nmatdn(ridx,i);
      }
      else if (i>ridx){
        Qtildeup(i-1) = Nmatup(i,ridx);
        Qtildedn(i-1) = Nmatdn(i,ridx);
        Rtildeup(i-1) = Nmatup(ridx,i);
        Rtildedn(i-1) = Nmatdn(ridx,i);
      }
    }
  
    const fptype paccup = 1.0/Stildeup;
    const fptype paccdn = 1.0/Stildedn;
    const fptype pacc = min(1.0, fptype(po)/p.K/fabs(paccup*paccdn)); //here prem^-1 occurs!, po and not po+1 appears because we are doing n+1 -> n
    
    #if DEBUG
    const fptype detninvup_np1 = Nmatup.inverse().determinant();
    const fptype detninvdn_np1 = Nmatdn.inverse().determinant();
    #endif
    
    const fptype r = rng_ptr->get_value();
    if(r<pacc){ //accept update
      update_accepted = 1;
      config_ptr->remove(ridx);
      //calculate P tilde matrix
      Eigen::MatrixXd Ptildeup = Eigen::MatrixXd::Zero(po-1, po-1);
      Eigen::MatrixXd Ptildedn = Eigen::MatrixXd::Zero(po-1, po-1);
      if(0==ridx){
        Ptildeup.noalias() = Nmatup.block(1,1,po-1,po-1);
        Ptildedn.noalias() = Nmatdn.block(1,1,po-1,po-1);
      }
      else if((po-1)==ridx){
        Ptildeup.noalias() = Nmatup.block(0,0,po-1,po-1);
        Ptildedn.noalias() = Nmatdn.block(0,0,po-1,po-1);
      }
      else{ //general case with four blocks
        //upper left block
        Ptildeup.block(0,0,ridx,ridx).noalias() = Nmatup.block(0,0,ridx,ridx);
        Ptildedn.block(0,0,ridx,ridx).noalias() = Nmatdn.block(0,0,ridx,ridx);
        //lower left block
        Ptildeup.block(ridx,0,po-ridx-1,ridx).noalias() = Nmatup.block(ridx+1,0,po-ridx-1,ridx);
        Ptildedn.block(ridx,0,po-ridx-1,ridx).noalias() = Nmatdn.block(ridx+1,0,po-ridx-1,ridx);
        //upper right block
        Ptildeup.block(0,ridx,ridx,po-ridx-1).noalias() = Nmatup.block(0,ridx+1,ridx,po-ridx-1);
        Ptildedn.block(0,ridx,ridx,po-ridx-1).noalias() = Nmatdn.block(0,ridx+1,ridx,po-ridx-1);
        //lower right block
        Ptildeup.block(ridx,ridx,po-ridx-1,po-ridx-1).noalias() = Nmatup.block(ridx+1,ridx+1,po-ridx-1,po-ridx-1);
        Ptildedn.block(ridx,ridx,po-ridx-1,po-ridx-1).noalias() = Nmatdn.block(ridx+1,ridx+1,po-ridx-1,po-ridx-1);
      } 
    
      //update N matrices
      Nmatup.resize(po-1, po-1);
      Nmatdn.resize(po-1, po-1);
      Nmatup.noalias() = Eigen::MatrixXd::Zero(po-1, po-1);
      Nmatdn.noalias() = Eigen::MatrixXd::Zero(po-1, po-1);
      Nmatup.noalias() = Ptildeup - Qtildeup*Rtildeup.transpose()/Stildeup;
      Nmatdn.noalias() = Ptildedn - Qtildedn*Rtildedn.transpose()/Stildedn;
      
      #if DEBUG
      const fptype detninvup_n = Nmatup.inverse().determinant();
      const fptype detninvdn_n = Nmatdn.inverse().determinant();
      const fptype pacc_explicit = min(1.0, fptype(po)/p.K*detninvup_n*detninvdn_n/detninvup_np1/detninvdn_np1);
      cout << pacc << " " << pacc_explicit << endl;
      #endif
    }
  }//end if perturbation order > 0
}

void CTAUXSolver::measure_gf(){
  
  const int po = config_ptr->get_perturbation_order();
  const fptype randomtau = this->p.beta*rng_ptr->get_value();
 
  Eigen::VectorXd egup = Eigen::VectorXd::Zero(po);
  Eigen::VectorXd egdn = Eigen::VectorXd::Zero(po);
  for(int i=0;i<po;i++){
    egup(i) = egamma( 1, config_ptr->get_auxspin(i)) - 1.0;
    egdn(i) = egamma(-1, config_ptr->get_auxspin(i)) - 1.0;
  }
  
  Eigen::MatrixXd Mup = Nmatup*egup.asDiagonal();
  Eigen::MatrixXd Mdn = Nmatdn*egdn.asDiagonal();
  
  Eigen::VectorXd gfup = Eigen::VectorXd::Zero(po);
  Eigen::VectorXd gfdn = Eigen::VectorXd::Zero(po);
  for(int i=0;i<po;i++){
    const fptype modtau = config_ptr->get_time(i) + randomtau;
    gfup(i) = wfup_ptr->get_interpolated_value(modtau);
    gfdn(i) = wfdn_ptr->get_interpolated_value(modtau);
  }
  Eigen::VectorXd Sup = Mup*gfup;
  Eigen::VectorXd Sdn = Mdn*gfdn;
  
  //conventional binning measurement
  for(int i=0;i<po;i++){
    //modify tau so that reference point is sampled
    const fptype modtau = config_ptr->get_time(i) + randomtau;
    fptype modtaurem = fmod(modtau, this->p.beta);
    const int prefactor = (modtau > this->p.beta) ? -1 : 1;
    //bin measured data
    const int binidx = floor(modtaurem/binwidth);
    gfupbins[binidx] += prefactor*Sup(i)/fptype(this->p.nsamplesmeasure);
    gfdnbins[binidx] += prefactor*Sdn(i)/fptype(this->p.nsamplesmeasure);
  }
  
  //Legendre coefficient measurement
  for(int i=0;i<this->p.nlegendre;i++){
    const fptype prefactor = sqrt(2.0*i+1.0)/fptype(this->p.nsamplesmeasure);
    Eigen::VectorXd lvec = Eigen::VectorXd::Zero(po);
    for(int j=0;j<po;j++){
      lvec(j) = outputgflegendreup_ptr->legendre_p(i, outputgflegendreup_ptr->x(config_ptr->get_time(j)));
    }
    gfuplegendre[i] += prefactor*lvec.dot(Sup);
    gfdnlegendre[i] += prefactor*lvec.dot(Sdn);
  }
}

void CTAUXSolver::construct_interacting_gf(){
  //this function constructs the interacting gf from the binned data
  for(int i=0;i<this->p.nbins;i++){
    const fptype tau = outputgfup_ptr->get_time(i);
    fptype valup = 0;
    fptype valdn = 0;
    for(int j=0;j<this->p.nbins;j++){
      //integration is approximated by rectangles, the bin width drops out in integration
      valup += wfup_ptr->get_interpolated_value(tau - binmids[j])*gfupbins[j];
      valdn += wfdn_ptr->get_interpolated_value(tau - binmids[j])*gfdnbins[j];
    }
    valup += wfup_ptr->get_interpolated_value(tau);
    valdn += wfdn_ptr->get_interpolated_value(tau);
    outputgfup_ptr->set_value(i, valup);
    outputgfdn_ptr->set_value(i, valdn);
  }
  
  //this function takes care of writing measured Legendre coefficients back to output classes
  for(int i=0;i<this->p.nlegendre;i++){
    outputgflegendreup_ptr->set_coefficient(i, gfuplegendre[i]);
    outputgflegendredn_ptr->set_coefficient(i, gfdnlegendre[i]);
  }
}

void CTAUXSolver::measure_perturbation_order(){
  
  average_po += config_ptr->get_perturbation_order()/fptype(this->p.nsamplesmeasure);
}

void CTAUXSolver::measure_acceptance_ratio(){
  
  acceptance_ratio += update_accepted/fptype(this->p.nsamplesmeasure);
}

#if DEBUG
void CTAUXSolver::calculate_Ninverse(){
  //this function calculates the inverse N matrices directly from the configuration of auxiliary spins and the input Weiss field to check the fast update formulas
  const int po = config_ptr->get_perturbation_order();
  
  Eigen::VectorXd egvecup = Eigen::VectorXd::Zero(po);
  Eigen::VectorXd egvecdn = Eigen::VectorXd::Zero(po);
  Eigen::MatrixXd gfmatup = Eigen::MatrixXd::Zero(po, po);
  Eigen::MatrixXd gfmatdn = Eigen::MatrixXd::Zero(po, po);
  for(int i=0;i<po;i++){
    egvecup(i) = egamma( 1, config_ptr->get_auxspin(i));
    egvecdn(i) = egamma(-1, config_ptr->get_auxspin(i));
    for(int j=0;j<po;j++){
      gfmatup(i,j) = wfup_ptr->get_interpolated_value(config_ptr->get_time(i) - config_ptr->get_time(j));
      gfmatdn(i,j) = wfdn_ptr->get_interpolated_value(config_ptr->get_time(i) - config_ptr->get_time(j));
    }
  }
  Eigen::MatrixXd Ninvup = egvecup.asDiagonal();
  Ninvup -= (egvecup - Eigen::VectorXd::Ones(po)).asDiagonal()*gfmatup;
  Eigen::MatrixXd Ninvdn = egvecdn.asDiagonal();
  Ninvdn -= (egvecdn - Eigen::VectorXd::Ones(po)).asDiagonal()*gfmatdn; 
  
  Eigen::MatrixXd nmatprodup = Ninvup*Nmatup - Eigen::MatrixXd::Identity(po,po);
  Eigen::MatrixXd nmatproddn = Ninvdn*Nmatdn - Eigen::MatrixXd::Identity(po,po);
  
  fptype sum = 0;
  for(int i=0;i<po;i++){
    for(int j=0;j<po;j++){
      sum += fabs(nmatprodup(i,j)) + fabs(nmatproddn(i,j));
    }
  }
  cout << sum << endl; //output should be zero if fast update algorithm works properly
}
#endif