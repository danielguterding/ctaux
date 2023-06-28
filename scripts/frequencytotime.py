#converts an imaginary frequency Green's function into a imaginary time Green's function
import os
import sys
import numpy as np

def load_gf(infilename):
  infilehandle = open(infilename, 'r')
  lines = infilehandle.readlines()[1:]
  infilehandle.close()
  nlines = len(lines)
  
  matsubarafreq = np.zeros((nlines), dtype=float)
  gf = np.zeros((nlines), dtype=complex)
  
  for i,l in enumerate(lines):
    splitline = l.strip().split()
    matsubarafreq[i] = float(splitline[0])
    gf[i] = complex(float(splitline[1]), float(splitline[2]))
    
  return matsubarafreq, gf

def transform_gf(matsubarafreq, gf, ntau):
  #calculate beta from first matsubara frequency w(n=0) = pi/beta*(2n+1) = pi/beta
  beta = np.pi/matsubarafreq[0]
  
  tau = np.linspace(0,beta,num=ntau,endpoint=True)
  newgf = np.zeros((ntau))
  
  for i,t in enumerate(tau):
    newgf[i] = 2.0/beta*np.sum(np.real(gf)*np.cos(t*matsubarafreq) + (np.imag(gf) + np.power(matsubarafreq, -1))*np.sin(t*matsubarafreq)) - 0.5

  return tau, newgf, beta

def write_gf(outfilename, tau, gf_time):
  outfilehandle = open(outfilename, 'w')
  outfilehandle.write('#tau in eV^-1 from 0 to beta; GF in eV^-1 real part, imag part is zero;\n')
  ntau = len(tau)
  for i in range(ntau):
    t = tau[i]
    gf = gf_time[i]
    outfilehandle.write('%1.14e % 1.14e\n' % (t, gf))
  outfilehandle.close()

def main():
  if(4 == len(sys.argv)):
    infilename = sys.argv[1]
    outfilename = sys.argv[2]
    ntau = int(float(sys.argv[3]))
    
    matsubarafreq, gf_frequency = load_gf(infilename)
    tau, gf_time, beta = transform_gf(matsubarafreq, gf_frequency, ntau)
    write_gf(outfilename, tau, gf_time)
    
    print("Beta is: %f" % beta)
    print("Data converted successfully.")
  else:
    print("Wrong number of input arguments. Please supply input Green's function, output file name and number of intended imaginary time points.")
  
main()
