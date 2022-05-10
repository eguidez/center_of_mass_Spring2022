# This code computes the xyz coordinates of the center of mass of two molecules
## Author: Emilie Guidez Date: 04/2022

INPUT: The .xyz file containing the coordinates of the molecules.
       The atomic_mass file containing the atomic mass in amu of the atoms in the molecules
OUTPUT: file named output_center_of_mass containing the x,y,z coordinates of the center of mass of each molecules

EXECUTION: To execute the code enter: python dimer_geometries.py
As prompted, enter the number of atoms in each molecule and the name of the .xyz file containing the coordinates.

TEST_INPUT: The water.xyz file contains the coordinates of a water dimer in Angstrom.
TEST_OUTPUT (water_dimer_output):
Molecule 1 center of mass coordinates:  -1.5192235062480295 -0.056766187419790295 0.0
Molecule 2 center of mass coordinates: 1.387526052619776 0.05717480447264767 0.0
