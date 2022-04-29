#!/usr/bin/env python
# coding: utf-8

# In[27]:


#This code will compute the coordinates of the center of mass of two molecules
#coordinates at various distances
import os
import numpy as np
#Input files
file_name = 'water_dimer.xyz'
file_path = os.path.abspath(file_name)
file = open(file_path,'r')
data = file.readlines()
file.close()
#print(data)

atom_file_name = 'atomic_mass'
atom_file_path = os.path.abspath(atom_file_name)
atom_file = open(atom_file_path,'r')
atom_mass_data = atom_file.readlines()
atom_file.close()
#print(atom_mass_data)

#Input data
atoms_molecule1 = 3
atoms_molecule2 =3
min_distance = 2.0  #Angstrom
max_distance = 3.0  #Angstrom
interval = 0.1 #A

# Assign data from the files
atom_names = np.empty([0])
x_coords = np.empty([0])
y_coords = np.empty([0])
z_coords = np.empty([0])
for linenum,line in enumerate(data):
    if (linenum == 0):
        total_atom_num = int(line)
    elif (linenum == 1):
        comment_line = line
    else:
        atom_data = line.split()
        atom_names = np.append(atom_names,atom_data[0])
        x_coords = np.append(x_coords,float(atom_data[1]))
        y_coords = np.append(y_coords,float(atom_data[2]))
        z_coords = np.append(z_coords,float(atom_data[3]))
#print(x_coords)
#print(atom_names)
#Get the mass of each atom
atom_mass = np.empty([0])
for i in range(0,len(atom_names)):
    for atom in atom_mass_data :
        atom_split = atom.split()
#        print('atom_split',atom_split)
        if atom_names[i] == atom_split[0]:
            atom_mass = np.append(atom_mass,float(atom_split[1]))
#print(atom_mass)

x_coords_1 = x_coords[0:atoms_molecule1]
y_coords_1 = y_coords[0:atoms_molecule1]
z_coords_1 = z_coords[0:atoms_molecule1]
atom_names_1 = atom_names[0:atoms_molecule1]
atom_mass_1 = atom_mass[0:atoms_molecule1]
x_coords_2 = x_coords[atoms_molecule1:total_atom_num]
y_coords_2 = y_coords[atoms_molecule1:total_atom_num]
z_coords_2 = z_coords[atoms_molecule1:total_atom_num]
atom_names_2 = atom_names[atoms_molecule1:total_atom_num]
atom_mass_2 = atom_mass[atoms_molecule1:total_atom_num]
#print(atom_names_1)
#print(atom_names_2)
#print(atom_mass_1)
#print(x_coords_1)
sum_mass_times_x_1 =0 
sum_mass_times_y_1 =0 
sum_mass_times_z_1 =0 
mass_sum_1 = 0
#Calculate the center of mass of each molecule
for i in range (0,atoms_molecule1):
    mass_times_x_1 = atom_mass_1[i] * x_coords_1[i]
    mass_times_y_1 = atom_mass_1[i] * y_coords_1[i]
    mass_times_z_1 = atom_mass_1[i] * z_coords_1[i]
    sum_mass_times_x_1 =sum_mass_times_x_1 + mass_times_x_1
    sum_mass_times_y_1 =sum_mass_times_y_1 + mass_times_y_1
    sum_mass_times_z_1 =sum_mass_times_z_1 + mass_times_z_1
    mass_sum_1 = mass_sum_1 + atom_mass_1[i]

com_x_molecule1 = sum_mass_times_x_1/mass_sum_1
com_y_molecule1 = sum_mass_times_y_1/mass_sum_1
com_z_molecule1 = sum_mass_times_z_1/mass_sum_1
print(com_x_molecule1,com_y_molecule1,com_z_molecule1)

sum_mass_times_x_2 =0 
sum_mass_times_y_2 =0 
sum_mass_times_z_2 =0 
mass_sum_2 = 0
for i in range (0,atoms_molecule2):
    mass_times_x_2 = atom_mass_2[i] * x_coords_2[i]
    mass_times_y_2 = atom_mass_2[i] * y_coords_2[i]
    mass_times_z_2 = atom_mass_2[i] * z_coords_2[i]
    sum_mass_times_x_2 =sum_mass_times_x_2 + mass_times_x_2
    sum_mass_times_y_2 =sum_mass_times_y_2 + mass_times_y_2
    sum_mass_times_z_2 =sum_mass_times_z_2 + mass_times_z_2
    mass_sum_2 = mass_sum_2 + atom_mass_2[i]

com_x_molecule2 = sum_mass_times_x_2/mass_sum_2
com_y_molecule2 = sum_mass_times_y_2/mass_sum_2
com_z_molecule2 = sum_mass_times_z_2/mass_sum_2
print(com_x_molecule2,com_y_molecule2,com_z_molecule2)

output_name = 'output_center_of_mass'
output = open(output_name,'w')
string1 = 'Molecule 1 center of mass coordinates: '+' '+ str(com_x_molecule1) +' '+ str(com_y_molecule1) +' '+ str(com_z_molecule1)+'\n'
string2 = 'Molecule 2 center of mass coordinates:'+' '+ str(com_x_molecule2) +' '+ str(com_y_molecule2) +' '+ str(com_z_molecule2)
output.write(string1)
output.write(string2)
output.close()


# In[ ]:




