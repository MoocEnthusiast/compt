"""
this file contains the main part of the program
combines the reading from file (of the parameters, conversions, setting up the grids, performing calculation and writing back to file
"""


## imports, declarations
from IOHANDLING import *
from CONSTANTS import *
import scipy
import numpy as np

## reading of the parameters
params = read_data('parameters.csv');
print params

Ee_0    = float(params["Ee_0"]);    #initial (central )energy
Ee_N    = int(params["Ee_N"]);      #number of (central ) energies for which im computing the spectrum
Ee_D    = float(params["Ee_D"]);    #distance between 2 succsesive energies

THETA01_L   = np.pi/180.*30.;#float(params["THETA01_L"]) #scattering angle
THETA02_L   = float(params["THETA02_L"]) #measurement angle
A0_L        = float(params["A0_0"])
LAMBDA_L    = float(params["LAMBDA_0"]) * 1e-9;
ORD         = int(params["NLCS_ORD"]);      #maximum order for nlcs
OMEGA_L     = 2*np.pi * CONST_C / LAMBDA_L;

## conversions and other stuff
#grid size for the spectrum, having 
Ep_gridsize         = 1000;
#grid size for the electron distribution
Ee_distr_gridsize   = 500;
# fwhm for the electorn distribution (actually sigma. to be changed)
sigma               = 0.01;
ORD = 2;

#for each of the (central) energies
for i in range(Ee_N):
    Ee_central = Ee_0 + i*Ee_D
    #now create a distribution of energies
    Ee_sigma            = sigma*Ee_central;
    #99.7% of distribution included between +/- 3*sigma  
    #this is the energies range for the distribution
    Ee_distr_energies   = np.linspace(Ee_central - 3*Ee_sigma, Ee_central + 3*Ee_sigma, Ee_distr_gridsize );
    #this is the (gaussian) distribution of energies
    Ee_distr            = nlcs_gaussian(Ee_distr_energies, Ee_central, Ee_sigma);
    
    #take the maximum energy and plug it in into nlcs
    
    Ep_cs_max = 0;
    for j in range(1,ORD+1):
        #compute compton scattering maximum energy in this distribution, which should be Ee_distr_energies[len(Ee_distr_energies)-1];      
        GAMMA_E = 1e+9*Ee_distr_energies[len(Ee_distr_energies)-1]/(CONST_ME*CONST_C**2*CONV_J2eV);
        
        Ep_cs = nlcs_max_scatter_energy_full(ORD,GAMMA_E, OMEGA_L,THETA01_L,THETA02_L,A0_L);
        
        if (Ep_cs > Ep_cs_max) :
            Ep_cs_max = Ep_cs;
    
        print Ep_cs_max*CONV_J2eV*1e-6;    
        #now u have the maximum energy for the grid for cross sections
        #grid cross section of the size energies and order
        cross_section_grid = np.zeros((ORD,Ep_gridsize));
        
        for j in range(Ee_distr_gridsize):
            
            
            Ej  = Ee_distr_energies[j];
            Rhoj= Ee_distr[j] 
            GAMMA_Ej = 1e+9*Ej/(CONST_ME*CONST_C**2*CONV_J2eV);
            for k in range(ORD):
                cross_section_grid[k][:] = Rhoj*nlcs_get_cross_section(k+1,Ep_cs_max,Rhoj,Ep_gridsize,GAMMA_Ej)+\
                                            cross_section_grid[k][:];
             
        Ep_grid = np.linspace(0,Ep_cs_max*CONV_J2eV*1e-6,Ep_gridsize)        
        plot(Ep_grid,cross_section_grid[0][:]);#+cross_section_grid[1][:]);
        
""" for each central_energy 
     create distribution
     take maximumm energy
     for each order of nlcs
         compute maximum
     declare grid according to maximum
     
     #done
     for each energy in distribution
         calculate spectrum 
         [aka for each order of the nlcs
             compute spectrum]
             
         add to grid     
"""























































