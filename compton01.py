from matplotlib  import pyplot;
import scipy
import numpy as np
from scipy.special import jn

from CONSTANTS import *

CONV_eV2J= 1.602e-19*1.0; 
CONV_J2eV= 1/CONV_eV2J;

LAMBDA_L    = 1050 * 1e-09;
OMEGA_L     = 2*np.pi * CONST_C / LAMBDA_L
GAMMA_E     = 99500;
BETA_E      = np.sqrt( (GAMMA_E**2.-1.) /GAMMA_E**2. )
THETA0_L    = np.pi/180.*17.;
THETA1_L    = 0.;
A0_L        = 0.4;
N_ORDER     = 1;



FINAL_ENERGY = 2* N_ORDER * OMEGA_L *CONST_HBAR * GAMMA_E**2 * (1 + BETA_E*np.cos(THETA0_L)) / \
               (\
               2*GAMMA_E**2* (1- BETA_E*np.cos(THETA1_L)) + \
               (2*N_ORDER*GAMMA_E* CONST_HBAR *OMEGA_L/(CONST_ME*CONST_C**2) + A0_L**2/(1+BETA_E*np.cos(THETA0_L))) * (1+np.cos(THETA0_L+THETA1_L)) \
               )

MAXIM_ENERGY = 2*N_ORDER *CONST_HBAR * OMEGA_L *GAMMA_E**2 * (1+np.cos(THETA0_L))/\
            (1+2*N_ORDER*CONST_HBAR*OMEGA_L*GAMMA_E/(CONST_ME*CONST_C**2)*(1+np.cos(THETA0_L)) +A0_L**2)
#energy in GeV or CONST_MEV
Ef = FINAL_ENERGY * CONV_J2eV * 1e-9
print "Maximum energy 01:",Ef
print "Maximum energy 02:",MAXIM_ENERGY*CONV_J2eV * 1e-9
theta = np.linspace(0,np.pi*2,720);

KN = (1+np.cos(theta)**2)/2*(1+GAMMA_E*(1-np.cos(theta))**2)**(-1) *\
    (GAMMA_E**2*(1-np.cos(theta))**2/((1+np.cos(theta)**2)*(1+GAMMA_E*(1-np.cos(theta))))+1)
#plt(KN);
N_ENERG = 1000
energy      = np.linspace(0,FINAL_ENERGY,N_ENERG);
cross_sect  = np.linspace(0,1, N_ENERG);

for i in range(0,N_ENERG):

    ENERG_E = GAMMA_E*CONST_ME*CONST_C**2;
    rho = 1;
    #initialize varibales    
    #u1
    u1 = 2.* (ENERG_E) * OMEGA_L * CONST_HBAR * (1. - np.cos(np.pi - THETA0_L)) / (CONST_ME**2 *CONST_C**4 * (1.+A0_L**2.))
    #u
    u  = energy[i] / (ENERG_E - energy[i]) 
    #un
    un = N_ORDER * u1 *1.
    #z
    z = 2. * A0_L / u1 * np.sqrt( u * (un - u) / ( 1. + A0_L**2  ) )
    
    CONST =  rho *CONST_QE**2.*CONST_ME**2.*CONST_C**4. / (16.*np.pi*ENERG_E) ;
    CONST = 1
    cross_sect[i]  = CONST * ( -4.* jn(N_ORDER,z)**2. + A0_L**2.*(2.+u**2./(1+u) )* \
    ( jn(N_ORDER-1,z)**2. + jn(N_ORDER+1,z)**2. - 2.*jn(N_ORDER,z)**2. ) )       
    
    #print cross_sect[i]; 
pyplot.plot(energy*CONV_J2eV * 1e-9,cross_sect)
#pyplot.yscale('log'    )










    