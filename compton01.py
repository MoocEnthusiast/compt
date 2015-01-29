from matplotlib.pyplot import plot as plt;
import scipy
import numpy as np
from scipy.special import jn

H   = 6.626e-34;
C   = 2.997e+08;
ME  = 9.109e-31;
QE  = 1.602e-19;
HBAR=H/(2*np.pi);

eV2J= 1.602e-19*1.0; 
J2eV= 1/eV2J;

LAMBDA_L    = 800 * 1e-09;
OMEGA_L     = 2*np.pi * C/ LAMBDA_L
GAMMA_E     = 3900;
BETA_E      = np.sqrt( (GAMMA_E**2.-1.) /GAMMA_E**2. )
THETA0_L    = np.pi/6;
THETA1_L    = 0;
A0_L        = 2;
N_ORDER     = 1;



FINAL_ENERGY = 2* N_ORDER * OMEGA_L * HBAR * GAMMA_E**2 * (1+np.cos(THETA0_L)) / \
               (\
               2*GAMMA_E**2* (1- BETA_E*np.cos(THETA1_L)) + (2*N_ORDER*GAMMA_E*HBAR*OMEGA_L/ME/C**2 + A0_L**2/(1+BETA_E*np.cos(THETA0_L))) * (1+np.cos(THETA0_L+THETA1_L)) \
               )

#energy in GeV or MeV
Ef = FINAL_ENERGY * J2eV * 1e-6 

theta = np.linspace(0,np.pi*2,720);

KN = (1+np.cos(theta)**2)/2*(1+GAMMA_E*(1-np.cos(theta))**2)**(-1) *\
    (GAMMA_E**2*(1-np.cos(theta))**2/((1+np.cos(theta)**2)*(1+GAMMA_E*(1-np.cos(theta))))+1)
#plt(KN);
N_ENERG = 1000
energy      = np.linspace(0,FINAL_ENERGY,N_ENERG);
cross_sect  = np.linspace(0,1, N_ENERG);

for i in range(0,N_ENERG):

    ENERG_E = GAMMA_E*ME*C**2;
    rho = 1;
    #initialize varibales    
    #u1
    u1 = 2.* (ENERG_E) * OMEGA_L * H * (1. - np.cos(2*np.pi - THETA0_L)) / (ME**2*C**4. * (1.+A0_L**2.))
    #u
    u  = energy[i] / (ENERG_E - energy[i]) 
    #un
    un = N_ORDER * u1 *1.
    #z
    z = 2. * A0_L / u1 * np.sqrt( u * (un - u) / ( 1. + A0_L**2  ) )
    
    CONST =  rho * QE**2.*ME**2.*C**4. / (16.*np.pi*ENERG_E) ;
    CONST = 1
    cross_sect[i]  = CONST * ( -4.* jn(N_ORDER,z)**2. + A0_L**2.*(2.+u**2./(1+u) )* \
    ( jn(N_ORDER-1,z)**2. + jn(N_ORDER+1,z)**2. - 2.*jn(N_ORDER,z)**2. ) )       
    
    print cross_sect[i]; 
plt(cross_sect)
    