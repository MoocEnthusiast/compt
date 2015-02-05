# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 11:45:33 2015

@author: calinh
"""
"""
Reads from a csv file the parameters for simulation
Values are comma separated


"""
def read_data(filename):
    paramsfile = open(filename,'r')
    param_list = [];
    param_dict = {};
    index = 0;
    for line in paramsfile:
        index = index + 1;
        
        l = line[:len(line)-1].split(',');   
        n = len(l);
        if index>4:
            n = len(l)
        for j in range(n):
            
            if (j%2==1):
                print float(l[j])
                param_list.append(float(l[j]))
                param_dict[l[j-1]]=float(l[j]);
                
        #print line
        #print param_dict
    return param_dict;
    
#param = get_params_from_file('parameters.csv');

#for x in param:
#    print x,"==>",param[x]

def write_data(Ee_central,Ee_fwhm,A0,order,energy_array,cross_section,other_parameters):
    
    # define filename
    # relevant quantities : electron energy, electron spread, a0, order of scattering
    # the rest will be written to file, this will appear in the name
    filename = "Ee="+str(Ee_central)+\
        "FWHM_Ee="+str(Ee_fwhm)+\
        "A0="+str(A0)+\
        "ORD="+str(order)+".txt"        
    data_file = open('filename','w')
    
    for x in other_parameters:
        data_file.write()
    for i in range(N_ENERG):
        #str_word=str(photon_energies[i]*1e+06)+" "+str(cross_section[i]/max(cross_section))+"\n"
   
       data_file.write(
        "{0:15.2f}  {1:.5f} \n".format(
            photon_energies[i]*1e+06,
            cross_section[i]/max(cross_section)
            ) 
        )
    data_file.close();
    
    
    