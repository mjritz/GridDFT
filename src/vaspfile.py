import numpy as np
import copy
import string
import sys


def readvasp_outcar(filename):
    forcexyz=[]
    xyz=[]
    final_energy=0
    opt_done=False
    coord_done=False
    f=open(filename,'r')
    vasp_outcar=f.readlines()
    i=np.size(vasp_outcar)-1

    #find the stationary point coordinates in the file
    while i > 0:
        if vasp_outcar[i].strip().startswith('POSITION'):
            #skip i ahead to coordinates
            if coord_done: break
            coord_start=copy.copy(i)+2
            coord_done=True
        i-=1

    i=np.size(vasp_outcar)-1
    while i > 0:
        if vasp_outcar[i].strip().startswith('free  energy   TOTEN'):
            final_energy=float(vasp_outcar[i].split()[4])
        i-=1

    #read coordinates
    i=coord_start
    while i <= np.size(vasp_outcar):
        if vasp_outcar[i].strip().startswith('--'):break
        xyz.append(vasp_outcar[i].split()[0:3])
        i+=1

    #read forces
    i=coord_start
    while i <= np.size(vasp_outcar):
        if vasp_outcar[i].strip().startswith('--'):break
        forcexyz.append(vasp_outcar[i].split()[3:6])
        i+=1
        #f.close()

    #read lattice vectors
    i=np.size(vasp_outcar)-1
    while i > 0:
        if vasp_outcar[i].strip().startswith('A1'):
           A1xx=str(((vasp_outcar[i].split()[3])))
           A1x=float(A1xx[:9])
           A1yy=str(vasp_outcar[i].split()[4])
           A1y=float(A1yy[:10])
           A1zz=str(vasp_outcar[i].split()[5])
           A1z=float(A1zz[:10])
           A2xx=str(vasp_outcar[i+1].split()[3])
           A2x=float(A2xx[:10])
           A2yy=str(vasp_outcar[i+1].split()[4])
           A2y=float(A2yy[:10])
           A2zz=str(vasp_outcar[i+1].split()[5])
           A2z=float(A2zz[:10])
           A3xx=str(vasp_outcar[i+2].split()[3])
           A3x=float(A3xx[:10])
           A3yy=str(vasp_outcar[i+2].split()[4])
           A3y=float(A3yy[:10])
           A3zz=str(vasp_outcar[i+2].split()[5])
           A3z=float(A3zz[:10])
        i-=1

    x_periodic = abs(np.array([A1x, A1y, A1z]))
    y_periodic =abs(np.array([A2x, A2y, A2z]))
    z_periodic = abs(np.array([A3x, A3y, A3z]))

    return np.array(xyz), np.array(forcexyz,float), y_periodic, x_periodic, z_periodic, final_energy

def readvasp_poscar(filename):
    f=open(filename,'r')
    vasp_poscar=f.readlines()
    j=np.size(vasp_poscar)-1
    xyz_old=[]
    x_vector=[]
    y_vector=[]
    z_vector=[]
    atnames=[]
    atnumbers=[]
    atom_name_array=[]


    while j > 0:
        if vasp_poscar[j].strip().startswith('Cartesian'):
                old_coord_start=copy.copy(j)+1
        j-=1

        j=np.size(vasp_poscar)-1
        while j > 0:
                if vasp_poscar[j].strip().startswith('Direct'):
                       old_coord_start=copy.copy(j)+1
                j-=1


        atnames.append(vasp_poscar[5].split())
        atnames_array=np.array(atnames).T
        atnumbers.append(vasp_poscar[6].split())
        atnumber_array=np.array(atnumbers,dtype=int).T
        x_vector.append(vasp_poscar[2].split()[0:3])
        x_array=(np.array(x_vector, dtype=float)).T
        y_vector.append(vasp_poscar[3].split()[0:3])
        y_array=np.array(y_vector,dtype=float).T
        z_vector.append(vasp_poscar[4].split()[0:3])
        z_array=np.array(z_vector, dtype=float).T


        for i in range(np.size(atnames)):
                for j in range(atnumber_array[i]):
                        atom_name_array.append(atnames_array[i])

    #read old coordinates
        i=old_coord_start
        while i <= np.size(vasp_poscar):
                if vasp_poscar[i].strip()=='':break
                xyz_old.append(vasp_poscar[i].split()[0:3])
                i+=1
        return atnames_array, atom_name_array,  np.array(xyz_old)

def makeinput(template_filename,filename,xyz,z_frozen):
    #read in template
    f_template=open(template_filename,'r')
    f=open(filename,'w')
    packmol_template=string.Template(f_template.read())

    #make formated xyz string
    coord_str=""
    for (i, (x,y,z)) in enumerate(xyz):
        if z >z_frozen:
            coord_str += "%15.5f%15.5f%15.5f %s %s %s\n" % (x, y, z, 'T', 'T', 'T')
        if z <z_frozen:
            coord_str += "%15.5f%15.5f%15.5f %s %s %s\n" % (x, y, z, 'F', 'F', 'F')

    input_str=packmol_template.substitute({'coordinates':coord_str})

    f.write(input_str)
    f.close()
    f_template.close()

#SPECIFIC TO BUILDING A SURFACE GRID. 

def ligand_vector(molecule_filename):
        f = open(molecule_filename, 'r')
        vasp_poscar=f.readlines()
        j=np.size(vasp_poscar)-1
        ligands=[]
        x_vector=[]
        y_vector=[]
        z_vector=[]
        atnames=[]
        atnumbers=[]

        #Determine line to start reading coordinates
        while j>0:
                if vasp_poscar[j].strip().startswith('Direct'):
                        molecule_atoms_start=copy.copy(j)+1
                j-=1

        atnames.append(vasp_poscar[5].split())
        atnumbers.append(vasp_poscar[6].split())
        x_vector.append(vasp_poscar[2].split()[0:3])
        x_array=(np.array(x_vector, dtype=float)).T
        y_vector.append(vasp_poscar[3].split()[0:3])
        y_array=np.array(y_vector,dtype=float).T
        z_vector.append(vasp_poscar[4].split()[0:3])
        z_array=np.array(z_vector, dtype=float).T

        #Read in coordinates for Precursor
        i=molecule_atoms_start
        while i<=np.size(vasp_poscar):
                if vasp_poscar[i].strip()=='':break
                ligands.append(vasp_poscar[i].split()[0:3])
                i+=1
        

        #Convert direct coordinates to cartesian using CONTCAR
        ligand_array_direct=np.array(ligands,dtype=float)
        ligand_array_cartesian=np.zeros_like(ligand_array_direct,dtype=float)

        for i in range(ligand_array_direct.shape[0]):
                ligand_array_cartesian[i,0]=(ligand_array_direct[i,0]*x_array[0,0])+(ligand_array_direct[i,0]*x_array[1,0])+(ligand_array_direct[i,0]*x_array[2,0])
                ligand_array_cartesian[i,1]=(ligand_array_direct[i,1]*y_array[0,0])+(ligand_array_direct[i,1]*y_array[1,0])+(ligand_array_direct[i,1]*y_array[2,0])
                ligand_array_cartesian[i,2]=(ligand_array_direct[i,2]*z_array[0,0])+(ligand_array_direct[i,2]*z_array[1,0])+(ligand_array_direct[i,2]*z_array[2,0])

        #Set up center molecule as first set of coordinates
        center_ligand=ligand_array_cartesian[0]
        ligand_vectors=np.zeros([ligand_array_cartesian.shape[0]-1, 3])

        #create array with ligand vectors from center molecule
        for i in range(ligand_array_cartesian.shape[0]-1):
               ligand_vectors[i]=ligand_array_cartesian[i+1]-ligand_array_cartesian[0]
        return ligand_vectors

def precursor_grid_build(surface_filename, grid_size,distance_from_surface):
        f = open(surface_filename, 'r')
        vasp_poscar=f.readlines()
        j=np.size(vasp_poscar)-1
        x_vector=[]
        y_vector=[]
        z_vector=[]
        surface=[]
        atnames=[]
        atnumbers=[]

        #Determine line to start reading coordinates
        while j>0:
                if vasp_poscar[j].strip().startswith('Direct'):
                        surface_atoms_start=copy.copy(j)+1
                j-=1


        atnames.append(vasp_poscar[5].split())
        atnumbers.append(vasp_poscar[6].split())
        x_vector.append(vasp_poscar[2].split()[0:3])
        x_array=(np.array(x_vector, dtype=float)).T
        y_vector.append(vasp_poscar[3].split()[0:3])
        y_array=np.array(y_vector,dtype=float).T
        z_vector.append(vasp_poscar[4].split()[0:3])
        z_array=np.array(z_vector, dtype=float).T
        
        #Read in coordinates for Surface
        i=surface_atoms_start
        while i<=np.size(vasp_poscar):
                if vasp_poscar[i].strip()=='':break
                surface.append(vasp_poscar[i].split()[0:3])
                i+=1

        #Convert direct coordinates to cartesian using CONTCAR
        surface_array_direct=np.array(surface,dtype=float)
        surface_array_cartesian=np.zeros_like(surface_array_direct,dtype=float)

        for i in range(surface_array_direct.shape[0]):
                surface_array_cartesian[i,0]=(surface_array_direct[i,0]*x_array[0,0])+(surface_array_direct[i,0]*x_array[1,0])+(surface_array_direct[i,0]*x_array[2,0])
                surface_array_cartesian[i,1]=(surface_array_direct[i,1]*y_array[0,0])+(surface_array_direct[i,1]*y_array[1,0])+(surface_array_direct[i,1]*y_array[2,0])
                surface_array_cartesian[i,2]=(surface_array_direct[i,2]*z_array[0,0])+(surface_array_direct[i,2]*z_array[1,0])+(surface_array_direct[i,2]*z_array[2,0])


        #Create Grid with N coordinates at set Z value
        Max_Z_coordinate=np.max(surface_array_cartesian[:,2])
        Z_Grid_Value=Max_Z_coordinate+distance_from_surface

        x = np.linspace(0,1,num=grid_size, endpoint=False)
        y = np.linspace(0,1,num=grid_size, endpoint=False)


        gridx,gridy=np.meshgrid(x,y)
        Grid_Array=np.zeros([grid_size**2,3],float)

        for i in range(Grid_Array.shape[0]):
                Grid_Array[i,2]=Z_Grid_Value

        i=0
        #for i in range(grid_array.shape[0]):
        for k in range(grid_size):
                for j in range(grid_size):
                        Grid_Array[i,0]=(gridx[j,k]*x_array[0,0])+(gridx[j,k]*x_array[1,0])+(gridx[j,k]*x_array[2,0])
                        i+=1
        i=0
        #for i in range(grid_array.shape[0]):
        for k in range(grid_size):
                for j in range(grid_size):
                        Grid_Array[i,1]=(gridy[j,k]*y_array[0,0])+(gridy[j,k]*y_array[1,0])+(gridy[j,k]*y_array[2,0])
                        i+=1

        return surface_array_cartesian,Grid_Array

def PES_geom_build(surface_array_cartesian, Grid_Array, ligand_vectors, precursor_atom_name, surface_atom_name):
        precursor=np.zeros([Grid_Array.shape[0], ligand_vectors.shape[0]+1, 4])
        poscar_atoms=np.concatenate([precursor_atom_name, surface_atom_name])
        a,b=np.unique(poscar_atoms,return_inverse=True)
        c=np.array(b)
        for i in range(Grid_Array.shape[0]):
                for j in range(ligand_vectors.shape[0]+1):
                        if j==0:
                                precursor[i,j,:] =[c[j],Grid_Array[i,0], Grid_Array[i,1], Grid_Array[i,2]]
                        else:
                                precursor[i,j,:]=[c[j],Grid_Array[i,0]+ligand_vectors[j-1,0],Grid_Array[i,1]+ligand_vectors[j-1,1],Grid_Array[i,2]+ligand_vectors[j-1,2]]

        #PA = poscar_atoms[0:ligand_vectors.shape[0]]
        #SA = poscar_atoms[ligand_vectors.shape[0]+1:c.size]
        SA_1=c[ligand_vectors.shape[0]+1:c.size]
        SA=np.transpose(SA_1)
        Surface_Coord_Name=np.column_stack((SA,surface_array_cartesian))

        #creates tuple with size of [grid_shape^2, surface + precursor atoms, 4] 
        PES_Geom = np.zeros([Grid_Array.shape[0], poscar_atoms.shape[0], 4])
        PES_Geometry=np.zeros([Grid_Array.shape[0], poscar_atoms.shape[0], 3])
        for i in range(Grid_Array.shape[0]):
            PES_Geom[i]=np.concatenate((precursor[i],Surface_Coord_Name),0)
            PES_Geom[i]=PES_Geom[i][PES_Geom[i,:,0].argsort()]
            PES_Geometry[i]=np.delete(PES_Geom[i], 0,axis=1)
        
        return PES_Geometry
