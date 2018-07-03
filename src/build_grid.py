import numpy as np
import vaspfile
import sys
import os
import subprocess


#Do not need to change. Need to save files "surface_CONTCAR", "surface_OUTCAR" and "precursor_OUTCAR" in /input folder. 
#Surface and precursor geometry files must be in direct coordinates
#The first atom in the precursor file has to be the center atom. 
poscar_file = 'POSCAR'
surface_contcar = '../inputs/silica_CONTCAR'
surface_outcar = '../inputs/silica_OUTCAR'
precursor_file='../inputs/precursor_CONTCAR'
#Change to filename that is appropriate for system
folder='Grid_VDW'
#Build a grid that has 3 by 3 sites for precursor. 
Grid_Axis = 4
Distance_from_surface = 2

#Creates a newfolder under ../ with folder name set above. 
newstep_folder='../Grid_VDW/'
subprocess.call('mkdir '+newstep_folder,shell=True)

#Converts the surface_CONTCAR file to cartesian coordinates and...
surface_array_cartesian, Grid_Array=vaspfile.precursor_grid_build(surface_contcar,Grid_Axis, Distance_from_surface)

#Reads in direct coordinates from precursor_CONTCAR, converts them to cartesian and builds an array of ligand vectors based on the center atom (first coordinate line). 
ligand_vectors=vaspfile.ligand_vector(precursor_file)


surface_atoms, surface_atom_name,xyz_old=vaspfile.readvasp_poscar(surface_contcar)
precursor_atoms, precursor_atom_name, xyz_old =vaspfile.readvasp_poscar(precursor_file)
xyz_react,force_react,y_periodic,x_periodic,z_periodic,final_energy=vaspfile.readvasp_outcar(surface_outcar)

poscar_atoms_repeat =np.concatenate([surface_atoms, precursor_atoms])
poscar_atoms=np.unique(poscar_atoms_repeat)


PES_Geometry = vaspfile.PES_geom_build(surface_array_cartesian, Grid_Array, ligand_vectors, precursor_atom_name, surface_atom_name)
print PES_Geometry

for i in range((Grid_Axis*Grid_Axis)):
    os.system('mkdir ../%s/0%i' %( folder, i))
    os.system('cp ../inputs/vasp_files/* ../%s/0%i/.' %(folder,i))
    print PES_Geometry[i]
    vaspfile.makeinput('../inputs/POSCAR_template',poscar_file, PES_Geometry[i],10)
    os.system('mv POSCAR ../%s/0%i/.' %(folder, i))
    os.chdir('../%s/0%i/' %(folder, i))
    os.system('bsub  -J \"surface-%i\" < bvasp' %(i))
    os.chdir('../../src/')


#cat /gpfs_share/santiso/mjritz/VASP/Potential/F/POTCAR /gpfs_share/santiso/mjritz/VASP/Potential/H/POTCAR /gpfs_share/santiso/mjritz/VASP/Potential/O/POTCAR /gpfs_share/santiso/mjritz/VASP/Potential/Si/POTCAR /gpfs_share/santiso/mjritz/VASP/Potential/W/POTCAR > POTCAR
