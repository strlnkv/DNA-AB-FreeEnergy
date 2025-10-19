# DNA-AB-FreeEnergy
Input data and scripts for OpenMM and LAMMPS used to obtain results in the paper "Cation-DNA outer sphere coordination in DNA polymorphism"

CCGGGCCCGG (10mers folder) or CCGGGCCCGGCCGGGCCCGG oligomer (20mers folder) in 80% v/v alcohol solutions or pure water

initialize.py - pre-equilibrate the system

gensmp.py - generate samples for umbrella sampling

8-24 - folders with the generated samples and scripts to extend trajectories

run.py - for 20mers in water and methanol only calculate trajectories of A and B/C forms, though A form is still obtained by pulling the system from B/C to A (thus the initial trajectory for A form most likely won't be in equilibrium, at least in methanol)

rst.pdb - input coordinates with DNA in model B form and solvent at random positions and orientations (packed in structure_and_coordinates.rar)

rst.psf - input structure files (packed in structure_and_coordinates.rar)

dna.txt - sample input file for LAMMPS used to calculate potential energy components for the whole system, as well as for ions in different pockets; ions corresponding to different pockets in the current frame should be assigned to separate groups (i* in the input file) and the same computes should be defined for each group

lmp_dat.txt - sample input structure file for LAMMPS, the coordinates and box sizes here should be replaced with the ones from the trajectory frame being analyzed (packed in structure_and_coordinates.rar)
