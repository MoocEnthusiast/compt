from matplotlib  import pyplot;
import scipy
import numpy as np
from scipy.special import jn

from CONSTANTS import *
from IOHANDLING import *


LAMBDA_L    = 800 * 1e-09;
OMEGA_L     = 2*np.pi * CONST_C / LAMBDA_L
GAMMA_E     = 99500;
THETA0_L    = np.pi/180.*30.;
THETA1_L    = 0.;
A0_L        = 1;
N_ORDER     = 1;
N_ENERG     = 1000

def nlcs_gamma2beta(GAMMA_E):
    return np.sqrt( (GAMMA_E**2 - 1.) / (GAMMA_E**2))

def nlcs_energy2gamma(GAMMA_E):
    return CONST_ME*CONST_C**2*GAMMA_E    
    
    
def nlcs_max_scatter_energy_full(N_ORDER,GAMMA_E,OMEGA_L,THETA0_L,THETA1_L,A0_L):
    BETA_E = nlcs_gamma2beta(GAMMA_E);    
    FINAL_ENERGY = 2* N_ORDER * OMEGA_L *CONST_HBAR * GAMMA_E**2 * (1 + BETA_E*np.cos(THETA0_L)) / \
               (\
               2*GAMMA_E**2* (1- BETA_E*np.cos(THETA1_L)) + \
                (2*N_ORDER*GAMMA_E* CONST_HBAR *OMEGA_L/(CONST_ME*CONST_C**2) \
                    + A0_L**2/(1+BETA_E*np.cos(THETA0_L))) * \
                (1+np.cos(THETA0_L+THETA1_L)) \
               )
    return FINAL_ENERGY  
    
def nlcs_max_scatter_energy_simplified(N_ORDER,GAMMA_E,OMEGA_L,THETA0_L,THETA1_L,A0_L):
    MAXIM_ENERGY = 2*N_ORDER * CONST_HBAR * OMEGA_L * GAMMA_E**2 * (1+np.cos(THETA0_L))/\
            (1+2*N_ORDER*CONST_HBAR*OMEGA_L*GAMMA_E/(CONST_ME*CONST_C**2)*(1+np.cos(THETA0_L)) +A0_L**2)
    return MAXIM_ENERGY
    
def nlcs_get_cross_section(N_ORDER,FINAL_ENERGY,rho_e,N_ENERG,GAMMA_E):
    
    cross_sect  = np.linspace(0,1, N_ENERG);    
    energy      = np.linspace(0,FINAL_ENERGY,N_ENERG);
    
    for i in range(0,N_ENERG):
        ENERG_E = GAMMA_E*CONST_ME*CONST_C**2;
        #rho_e = 1;
        #initialize varibales    
        #u1
        u1 = 2.* (ENERG_E) * OMEGA_L * CONST_HBAR * (1. - np.cos(np.pi - THETA0_L)) / (CONST_ME**2 *CONST_C**4 * (1.+A0_L**2.))
        #u
        u  = energy[i] / (ENERG_E - energy[i]) 
        #un
        un = N_ORDER * u1 *1.
        #z
        z = 2. * A0_L / u1 * np.sqrt( u * (un - u) / ( 1. + A0_L**2  ) )
        
        CONST =  rho_e*CONST_QE**2.*CONST_ME**2.*CONST_C**4. / (16.*np.pi*ENERG_E) ;
        CONST = 1
        cross_sect[i]  = CONST * ( -4.* jn(N_ORDER,z)**2. + A0_L**2.*(2.+u**2./(1+u) )* \
            ( jn(N_ORDER-1,z)**2. + jn(N_ORDER+1,z)**2. - 2.*jn(N_ORDER,z)**2. ) )       
    return cross_sect;

def nlcs_gaussian(x, mean, sig):
    return np.exp( - (x - mean)**2. / (2. * sig**2) )/ (sig* np.sqrt(2*np.pi))
    
##########################################
##########################################
#intermidiate calculation
##########################################
##########################################
    
FINAL_ENERGY = nlcs_max_scatter_energy_full(N_ORDER,GAMMA_E,OMEGA_L,THETA0_L,THETA1_L,A0_L);
MAXIM_ENERGY = nlcs_max_scatter_energy_simplified(N_ORDER,GAMMA_E,OMEGA_L,THETA0_L,THETA1_L,A0_L);

#energy in GeV or CONST_MEV
Ef = FINAL_ENERGY * CONV_J2eV * 1e-9
print "Maximum energy 01:",Ef
print "Maximum energy 02:",MAXIM_ENERGY*CONV_J2eV * 1e-9

# klein nishina cross section 
#theta = np.linspace(0,np.pi*2,720);
#KN = (1+np.cos(theta)**2)/2*(1+GAMMA_E*(1-np.cos(theta))**2)**(-1) *\
#    (GAMMA_E**2*(1-np.cos(theta))**2/((1+np.cos(theta)**2)*(1+GAMMA_E*(1-np.cos(theta))))+1)
#plt(KN);




#energy      = np.linspace(0,FINAL_ENERGY,N_ENERG);
##########################################
##########################################
##########################################
##########################################



energy_mean     = 2.500e+09;
FWHM            = 0.1;
energy_sig      = FWHM*energy_mean/(2.*np.sqrt(2*np.log(2)));



energies        = 500;
energy_elec     = np.linspace(energy_mean-4*energy_sig,energy_mean+4*energy_sig,energies)
distr_elec      = gaussian(energy_elec,energy_mean,energy_sig)
distr_elec      = distr_elec/(max(distr_elec))

#plot(energy_elec,distr_elec)
GAMMA_MAX = energy_elec[energies-1]/(CONST_ME*CONST_C**2 * CONV_J2eV)
ABS_MAX = nlcs_max_scatter_energy_full(N_ORDER,GAMMA_MAX,OMEGA_L,THETA0_L,THETA1_L,A0_L);

cross_sections = np.zeros((energies,N_ENERG))
"""
for j in range(distr_elec.size):
    RHO_E = distr_elec[j];
    
    GAMMA_E = energy_elec[j]/(CONST_ME*CONST_C**2 * CONV_J2eV)
    FINAL_ENERGY = nlcs_max_scatter_energy_full(N_ORDER,GAMMA_E,OMEGA_L,THETA0_L,THETA1_L,A0_L);
    NR_ELEM = int((FINAL_ENERGY/ABS_MAX) * N_ENERG)
    #print "Energy for step ",j," :",FINAL_ENERGY * CONV_J2eV * 1e-6
    cross_sections[j][:NR_ELEM] = RHO_E*nlcs_get_cross_section(N_ORDER,FINAL_ENERGY,RHO_E,NR_ELEM,GAMMA_E);
    #plot(test_cs)

final_cs = np.linspace(0,1,N_ENERG)*0.;
for j in range(distr_elec.size):
    final_cs[:]=final_cs[:]+cross_sections[j][:];
    
photon_energies = np.linspace(0,ABS_MAX*CONV_J2eV * 1e-6,N_ENERG)
plot(photon_energies,final_cs/max(final_cs))

filename = "A0="+str(A0_L)+\
            "|lambda="+str(LAMBDA_L*1e+9)+\
            "nm|Ee="+str(energy_mean*1e-9)+\
            "GeV|FWHM="+str(FWHM)+\
            "|ORD="+str(N_ORDER)+".txt"
            
data_file = open(filename,'w')

for i in range(N_ENERG):
    str_word=str(photon_energies[i]*1e+06)+" "+str(final_cs[i]/max(final_cs))+"\n"
   
    data_file.write(
    "{0:15.2f}  {1:.5f} \n".format(
            photon_energies[i]*1e+06,
            final_cs[i]/max(final_cs)
            ) 
        )

data_file.close();
    """
#pyplot.plot(energy*CONV_J2eV * 1e-9,test_cs)
#pyplot.yscale('log'    )










    