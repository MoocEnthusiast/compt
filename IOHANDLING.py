# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 11:45:33 2015

@author: calinh
"""
"""
Reads from a csv file the parameters for simulation
Values are comma separated


"""
import os;
def read_data(filename):
    paramsfile = open(filename,'r')
    param_dict = {};
    index = 0;
    
    for line in paramsfile:
        index = index + 1;
        
        l = line[:len(line)-1].split(',');   
        n = len(l);
        
        for j in range(n):
            if (j%2==1):
                param_dict[l[j-1]]=float(l[j]);
                
    return param_dict;
    
#param = get_params_from_file('parameters.csv');

#for x in param:
#    print x,"==>",param[x]

def write_data(Ee_central,Ee_fwhm,A0,order,photon_energies,cross_section,other_parameters):
    
    # define filename
    # relevant quantities : electron energy, electron spread, a0, order of scattering
    # the rest will be written to file, this will appear in the name
    if (os.path.isdir("./data")==False):
        os.mkdir("data");
        
    filename = "Ee="+str(Ee_central)+\
        "FWHM_Ee="+str(Ee_fwhm)+"%"+\
        "A0="+str(A0)+\
        "ORD="+str(order)+".txt"    
        
    directory = "./data"    
    data_file = open(directory+"/"+filename,'w')
    
    for x in sorted(other_parameters):
        data_file.write("{0} = {1} ;\n".format(x,other_parameters[x]))
        
    for i in range(len(photon_energies)):
        #str_word=str(photon_energies[i]*1e+06)+" "+str(cross_section[i]/max(cross_section))+"\n"
   
       data_file.write(
        "{0:15.2f}  {1:.5f} \n".format(
            photon_energies[i],
            cross_section[i]
            ) 
        )
    data_file.close();
    
    
    