#converts a density of states file into imaginary frequency Green's functions
import sys
import os
import numpy as np

def main():
  if(5 == len(sys.argv)):
    infilename = sys.argv[1]
    outfilename = sys.argv[2]
    nmatsubara = int(sys.argv[3])
    beta = float(sys.argv[4])
    
    infilehandle = open(infilename, 'r')
    lines = infilehandle.readlines()[1:] #discard comment line
    infilehandle.close()
    
    norbitals = len(lines[0].strip().split())-1
    ndosentries = len(lines)
    
    energies = np.zeros((ndosentries), dtype=float)
    dos = np.zeros((norbitals, ndosentries), dtype=float)
    
    for i,l in enumerate(lines):
      splitline = l.strip().split()
      energies[i] = splitline[0]
      for j,d in enumerate(splitline[1:]):
	dos[j, i] = d
   
    matsubarafrequencies = np.array([(2.*n+1.)*np.pi/beta for n in range(nmatsubara)])
    orbitalgf = np.zeros((norbitals, nmatsubara), dtype=complex)
    for i,orbitaldos in enumerate(dos):
      for j,wn in enumerate(matsubarafrequencies):
	complexintegrand = orbitaldos/(1j*wn - energies)
	realintegrand = np.real(complexintegrand)
	imagintegrand = np.imag(complexintegrand)
	orbitalgf[i, j] = np.trapz(y=realintegrand, x=energies) + 1j*np.trapz(y=imagintegrand, x=energies)
    
    for i,orbitalgf in enumerate(orbitalgf):
      outfilenamethisorb = outfilename + '_%0i' % i
      outfilehandle = open(outfilenamethisorb, 'w')
      outfilehandle.write('#matsubara frequencies in eV, real part of Greens function, imag part of Greens function\n')
      for wn,gf in zip(matsubarafrequencies, orbitalgf):
	outfilehandle.write('%1.14f % 1.14f % 1.14f\n' % (wn, gf.real, gf.imag))
      outfilehandle.close() 
  else:
    print 'Wrong number of input arguments. Please supply input DOS file, output file name, desired number of Matsubara frequencies and desired inverse Temperature beta in inverse eV.'
  
main()