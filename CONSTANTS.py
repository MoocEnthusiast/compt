from matplotlib.pyplot import plot
import scipy
import numpy as np


CONST_KB    = 1.38064881e-23
CONST_H     = 6.62606957e-34;
CONST_C     = 2.99792458e+08; 
CONST_ME    = 9.10938291e-31;
CONST_MP    = 1.67262178e-27;
CONST_QE    = 1.60217657e-19;
CONST_EPS   = 8.85418781e-12;
CONST_HBAR  = CONST_H/(2*np.pi);

CONV_eV2J   = CONST_QE * 1.;
CONV_J2eV   = 1./ CONV_eV2J;

 