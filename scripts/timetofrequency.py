#converts an imaginary time Green's function to an imaginary frequency Green's function
import os
import sys
import numpy as np

def load_gf(infilename):
  infilehandle = open(infilename, 'r')
  lines = infilehandle.readlines()[1:]
  infilehandle.close()
  nlines = len(lines)
  
  tau = np.zeros((nlines), dtype=float)
  gf = np.zeros((nlines), dtype=float)
  
  for i,l in enumerate(lines):
    splitline = l.strip().split()
    tau[i] = float(splitline[0])
    gf[i] = float(splitline[1])
    
  return tau, gf

def transform_gf(taupoints, gf, nfreq):
  beta = taupoints[-1]
  matsubarafreq = np.pi/beta*(2*np.array(range(nfreq))+1)
  newgf = np.zeros((nfreq), dtype=complex)
  
  for i in range(nfreq):
    newgf[i] = np.trapz(np.exp(1j*matsubarafreq[i]*taupoints)*gf,x=taupoints)
  
  return matsubarafreq, newgf, beta

def write_gf(outfilename, tau, gf_frequency):
  outfilehandle = open(outfilename, 'w')
  outfilehandle.write('#matsubara frequencies in eV, real part of Greens function, imag part of Greens function\n')
  ntau = len(tau)
  for i in range(ntau):
    t = tau[i]
    gf = gf_frequency[i]
    outfilehandle.write('%1.14e % 1.14e % 1.14e\n' % (t, gf.real, gf.imag))
  outfilehandle.close()

def main():
  if(4 == len(sys.argv)):
    infilename = sys.argv[1]
    outfilename = sys.argv[2]
    nfreq = int(sys.argv[3])
    
    taupoints, gf_time = load_gf(infilename)
    matsubarafreq, gf_frequency, beta = transform_gf(taupoints, gf_time, nfreq)
    write_gf(outfilename, matsubarafreq, gf_frequency)
    
    print("Beta is: %f" % beta)
    print("Data converted successfully.")
  else:
    print("Wrong number of input arguments. Please supply input Green's function, output file name and number of intended imaginary frequencies.")
  
main()