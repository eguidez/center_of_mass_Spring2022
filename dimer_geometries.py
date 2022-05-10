#!/usr/bin/env python
# coding: utf-8

# In[27]:


#This code will compute the coordinates of the center of mass of two molecules
#coordinates at various distances
import os
import numpy as np
#Input files
def read_file(file_name):
    """
    Reads file_name text file and returns the content in a data list
    INPUT: file_name
    OUTPUT: data
    """
    file_path = os.path.abspath(file_name)
    file = open(file_path,'r')
    data = file.readlines()
    file.close()
    return data

def assign_data(data):
    """
    Assign x,y,z coordinates and atom names from data to their respective array.
    INPUT: data array generated from reading an xyz file
    OUTPUT: x_coords, y_coords, z_coords, atom_names arrays
    """
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
# Generate atom_mass array here
    atom_mass = atom_mass_assign(atom_names)
    return (x_coords, y_coords, z_coords, atom_names,atom_mass)

def atom_mass_assign(atom_names):
    """
    Assign atomic mass in amu to each atom in atom_names array
    INPUT: Array atom_names containing the atomic symbols of the molecules
    OUTPUT: Array atom_mass containing the mass of each atom in the molecules
    """
    atom_mass = np.empty([0])
    for i in range(0,len(atom_names)):
        for atom in atom_mass_data :
            atom_split = atom.split()
#        print('atom_split',atom_split)
            if atom_names[i] == atom_split[0]:
                atom_mass = np.append(atom_mass,float(atom_split[1]))
    return (atom_mass)

def data_per_molecule(array,atom_number_array,molecule_number):
    """
    Slices array to only get the data for one molecule
    INPUT: array to slice (array), array containing number of atoms in
    each molecules (atom_number_array) and molecule number
    OUTPUT: array_for_one_molecule
    """
    if molecule_number == 0:
        start = 0
        end = int(atom_number_array[molecule_number])
    else :
        start = int(atom_number_array[molecule_number])
        end =  int(atom_number_array[0]) + int(atom_number_array[1])
    array_for_one_molecule = array[start:end]
    return (array_for_one_molecule)

def center_of_mass(x_coords,y_coords,z_coords,atom_mass):
    """
    Calculates center of mass coordinates of a molecule
    INPUT: x_coordinates (x_coords),y_coordinates(y_coord),
        z_coordinate(z_coord),atom masses (atom_mass)
    OUTPUT: com_x,com_y,com_z
    """
    sum_mass_times_x =0
    sum_mass_times_y =0
    sum_mass_times_z =0
    mass_sum = 0
#Calculate the center of mass of the molecule
    for i in range (0,len(atom_mass)):
        mass_times_x = atom_mass[i] * x_coords[i]
        mass_times_y = atom_mass[i] * y_coords[i]
        mass_times_z = atom_mass[i] * z_coords[i]
        sum_mass_times_x =sum_mass_times_x + mass_times_x
        sum_mass_times_y =sum_mass_times_y + mass_times_y
        sum_mass_times_z =sum_mass_times_z + mass_times_z
        mass_sum = mass_sum + atom_mass[i]

    com_x = sum_mass_times_x/mass_sum
    com_y = sum_mass_times_y/mass_sum
    com_z = sum_mass_times_z/mass_sum
    return(com_x,com_y,com_z)

def write_output(output_name,com_coords,molecule_number):
    """
    Write the center of masss coordinates of the molecule to output_name
    INPUT: output file name (output_name) and xyz coords of the center of mass
    OUTPUT: output file(output_name) containing xyz coords of the center of mass
    """
    output = open(output_name,'a')
    string = 'Molecule ' + str(molecule_number+1) + \
    ' center of mass coordinates: '+' '+ str(com_coords[0]) +' ' \
    + str(com_coords[1]) +' '+ str(com_coords[2])+'\n'
    output.write(string)
    output.close()

if __name__ == '__main__':
#Input data from the user
    atoms_molecule1 = input("Enter the number of atoms in molecule 1: ")
    atoms_molecule2 = input("Enter the number of atoms in molecule 2: ")
    file_name = input("Enter the name of the .xyz file containing the atom coordinates: ")
    output_name = input("Enter the name of the file containing the output ")

#Array containing the number of atoms in each molecule
    atom_number_array = np.array([atoms_molecule1, atoms_molecule2])

#Now extract data from the files
    xyz_data = read_file(file_name)
    atom_mass_data = read_file('atomic_mass')
    all_x, all_y, all_z, all_atom_names,all_atom_mass = assign_data(xyz_data)
    for molecule in range (0,len(atom_number_array)):
        x_coords = data_per_molecule(all_x,atom_number_array,molecule)
        y_coords = data_per_molecule(all_y,atom_number_array,molecule)
        z_coords = data_per_molecule(all_z,atom_number_array,molecule)
        atom_names = data_per_molecule(all_atom_names,atom_number_array,molecule)
        atom_mass = data_per_molecule(all_atom_mass,atom_number_array,molecule)
        com_coordinates = center_of_mass(x_coords,y_coords,z_coords,atom_mass)
        write_output(output_name,com_coordinates,molecule)
    print('The program terinated normally')
