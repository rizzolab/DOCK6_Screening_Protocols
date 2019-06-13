import math, sys
import os.path
import cmath
from math import sqrt

#################################################################################################################
#################################################################################################################
# data structure to store information about each residue with the docked ligand.
class Mol:
    def __init__(self,name,atom_list,bond_list,residue_list):
        self.name       = str(name)
        self.atom_list  = atom_list
        self.bond_list  = bond_list
	self.residue_list = residue_list

class atom:
    def __init__(self,X,Y,Z,Q,type,name,num,resnum,resname):
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
        self.Q = float(Q)
        self.heavy_atom = False
        self.type = type
        self.name = name
        self.num  = int(num)
	self.resnum  = int(resnum)
	self.resname = resname
class bond:
     def __init__(self,a1_num,a2_num,num,type):
        self.a1_num = int(a1_num)
        self.a2_num = int(a2_num)
        self.num = int(num)
        self.type = type
class residue:
     def __init__(self,atom_list,resnum,resname):
	self.atom_list = atom_list
	self.resnum  = int(resnum)
        self.resname = resname

#################################################################################################################
#################################################################################################################
def read_Mol2_file(file):
    # reads in data from multi-Mol2 file.

    file1 = open(file,'r')
    lines  =  file1.readlines()
    file1.close()

    atom_list = []
    bond_list = []
    residue_list = {}
    mol_list = []

    flag_atom    = False
    flag_bond    = False
    flag_substr  = False
    flag_mol     = False
    flag_getName = False

    i = 0  # i is the num of molecules read so far
    for line in lines:
         linesplit = line.split() #split on white space
         if (len(linesplit) == 1):
            if(linesplit[0] == "@<TRIPOS>MOLECULE"):
               i = i + 1
               #print "READING IN MOL #" + str(i)
               #print "read in molecule info:"
               line_num = 0
               flag_mol = True
               flag_atom = False
               flag_bond = False
               flag_substr = False

            if(linesplit[0] == "@<TRIPOS>ATOM"):
               #print "read in atom info:"
               flag_atom = True
               flag_bond = False
               flag_substr = False
               flag_mol = False

            if(linesplit[0] == "@<TRIPOS>BOND"):
               #print "read in bond info:"
               flag_bond = True
               flag_substr = False
               flag_mol = False
               flag_atom = False

            if(linesplit[0] == "@<TRIPOS>SUBSTRUCTURE"):
               #print "read in substructure info:"
               flag_substr = True
               flag_mol = False
               flag_atom = False
               flag_bond = False
         if (flag_mol and (not flag_getName) and len(linesplit)==1 ):
             if (line_num == 1):
                line_num = 0
                Name = linesplit[0]
                flag_getName = True
             line_num = line_num + 1

         if ((len(linesplit) >= 9 )and (flag_atom)):
             atom_num  = linesplit[0]
             atom_name = linesplit[1]
             X         = linesplit[2]
             Y         = linesplit[3]
             Z         = linesplit[4]
             atom_type = linesplit[5]
             res_num   = int(linesplit[6])
             res_name  = linesplit[7]
             Q         = linesplit[8]
             temp_atom = atom(X,Y,Z,Q,atom_type,atom_name,atom_num,res_num,res_name)
             atom_list.append(temp_atom)
	     if residue_list.has_key(res_num):
     		 residue_list[res_num].append(temp_atom)
    	     else:
                 residue_list[res_num] = [temp_atom]

         elif (len(linesplit) == 4 and flag_bond):
             bond_num  = linesplit[0]
             a1_num    = linesplit[1]
             a2_num    = linesplit[2]
             bond_type = linesplit[3]
             temp_bond = bond(a1_num,a2_num,bond_num,bond_type)
             bond_list.append(temp_bond)

         elif (flag_substr):
                 ID_heavy_atoms(atom_list)
                 data = Mol(Name,atom_list,bond_list,residue_list)
                 mol_list.append(data)
                 flag_getName = False
                 flag_substr = False
                 atom_list = [];bond_list = []

    return mol_list
#################################################################################################################
# Does not work with grid
#def write_mol2(molecule,filename):
#
#	outmol2 = open(filename,'w')
#	outmol2.write("@<TRIPOS>MOLECULE\n")      #start the MOLECULE RTI (Record Type Indicator)
#	outmol2.write(molecule.name+'\n')         #print MOL2FILE name of the molecule
#	outmol2.write(" %d %d %d 0 0\n" % (len(molecule.atom_list), 
#		len(molecule.bond_list), len(molecule.residue_list.keys()))) 
#	# For now, the number of residues is hard-coded to 1. To be fixed.
#	outmol2.write("SMALL\n") 		  #mol_type
#	outmol2.write("USER_CHARGES\n") 	  #charge_type
#
#	outmol2.write("\n@<TRIPOS>ATOM\n")      #start the ATOM RTI (Record Type Indicator)
#	for j in range(0,len(molecule.atom_list)):
#        	outmol2.write("%6d %-4s %9.4f %9.4f %9.4f %-5s %4s %6s %9.4f\n" % 
#		(j+1, molecule.atom_list[j].name, molecule.atom_list[j].X, molecule.atom_list[j].Y, 
#		molecule.atom_list[j].Z, molecule.atom_list[j].type, molecule.atom_list[j].resnum, 
#		molecule.atom_list[j].resname, molecule.atom_list[j].Q))
#
#	outmol2.write("@<TRIPOS>BOND\n")
#	for m in range(0,len(molecule.bond_list)):
#        	outmol2.write("%7d %5d %-5d %s\n" % (molecule.bond_list[m].num, 
#		molecule.bond_list[m].a1_num, molecule.bond_list[m].a2_num, molecule.bond_list[m].type))
#
#	outmol2.write("@<TRIPOS>SUBSTRUCTURE\n")
#	for resnum in molecule.residue_list.keys():
#		outmol2.write("%7d %8s %5d RESIDUE 1 A\n" % (resnum, 
#		molecule.residue_list[resnum][0].resname, # residue name 
#		molecule.residue_list[resnum][0].num ))   # atom num of first atom in this residue
#	outmol2.close()
#    	return
#################################################################################################################
#################################################################################################################
def write_mol2(molecule,filename):

        outmol2 = open(filename,'w')
        outmol2.write("@<TRIPOS>MOLECULE\n")      #start the MOLECULE RTI (Record Type Indicator)
        outmol2.write(molecule.name+'\n')         #print MOL2FILE name of the molecule
        outmol2.write(" %d %d %d 0 0\n" % (len(molecule.atom_list),
                len(molecule.bond_list), len(molecule.residue_list.keys())))
        # For now, the number of residues is hard-coded to 1. To be fixed.
        outmol2.write("SMALL\n")                  #mol_type
        outmol2.write("USER_CHARGES\n")           #charge_type

        #outmol2.write("\n@<TRIPOS>ATOM\n")      #start the ATOM RTI (Record Type Indicator)
        outmol2.write("@<TRIPOS>ATOM\n")      #start the ATOM RTI (Record Type Indicator)
        for j in range(0,len(molecule.atom_list)):
                outmol2.write("%-5d %-5s %9.4f %9.4f %9.4f %-5s %4s %-6s %8.4f\n" %
                (j+1, molecule.atom_list[j].name, molecule.atom_list[j].X, molecule.atom_list[j].Y,
                molecule.atom_list[j].Z, molecule.atom_list[j].type, molecule.atom_list[j].resnum,
                molecule.atom_list[j].resname, molecule.atom_list[j].Q))

        outmol2.write("@<TRIPOS>BOND\n")
        for m in range(0,len(molecule.bond_list)):
                outmol2.write("%-7d %5d %-5d %s\n" % (molecule.bond_list[m].num,
                molecule.bond_list[m].a1_num, molecule.bond_list[m].a2_num, molecule.bond_list[m].type))

        outmol2.write("@<TRIPOS>SUBSTRUCTURE\n")
        for resnum in molecule.residue_list.keys():
                outmol2.write("%-7d %8s %5d RESIDUE 1 A %3s 1\n" % (resnum,
                molecule.residue_list[resnum][0].resname, # residue name
                molecule.residue_list[resnum][0].num,   # atom num of first atom in this residue
                molecule.residue_list[resnum][0].resname[0:3] )) # residue
        outmol2.close()
        return
#################################################################################################################
def get_pdbcode_list(filename):
    systems_list = open(file,'r')
    lines  =  systems_list.readlines()
    return lines	  
#################################################################################################################
def ID_heavy_atoms(atom_list):
    for i in range(len(atom_list)):
        if (atom_list[i].type[0] != 'H'):
            atom_list[i].heavy_atom = True
    return atom_list
#################################################################################################################
#################################################################################################################
def distance2(vector1,vector2):
        distance2 = 0
        distance2 += (vector1.X-vector2.X)**2
        distance2 += (vector1.Y-vector2.Y)**2
        distance2 += (vector1.Z-vector2.Z)**2
        return distance2
#################################################################################################################
#################################################################################################################
def norm(vector1):
        norm = 0
        for i in range(len(vector1)):
                norm += (vector1[i])*(vector1[i])
        return sqrt(norm)
#################################################################################################################
#################################################################################################################
# Make sure the heavy atoms are being declared as heavy
# i.e call ID_heavy atoms function
def heavy_atom_RMSD(ref,pose):
    if (len(ref.atom_list) != len(pose.atom_list)):
       return -1 # when atom numbers do not agree
    sum = 0.0
    num_hvy_atoms = 0
    for i in range(len(ref.atom_list)):
        if (ref.atom_list[i].heavy_atom and pose.atom_list[i].heavy_atom):
           sum += distance2(ref.atom_list[i],pose.atom_list[i])
           num_hvy_atoms+=1
    return sqrt(sum/num_hvy_atoms)

#################################################################################################################
def formal_charge(molecule):
        total = 0
        for i in range(len(molecule.atom_list)):
                total += molecule.atom_list[i].Q
        return total
#################################################################################################################
def centre_of_mass(molecule):
        # Dictionary of atomic weights of elements
        atom_mass = {'O':15.9994 ,'N':14.00674 ,'C':12.011 ,'F':18.9984032 ,'Cl':35.4527 ,'Br':79.904
        ,'I':126.90447 ,'H':1.00794 ,'B':10.811 ,'S':32.066 ,'P':30.973762 ,'Li':6.941 ,'Na':22.98968
        ,'Mg':24.3050 ,'Al':26.981539 ,'Si':28.0855 ,'K':39.0983 ,'Ca':40.078 ,'Cr':51.9961 ,'Mn':54.93805
        ,'Fe':55.847 ,'Co':58.93320 ,'Cu':63.546 ,'Zn':65.39 ,'Se':78.96 ,'Mo':95.94 ,'Sn':118.710 ,'LP':0.0 }

        cmass = [0,0,0]
        centroid = [0,0,0]
        molecular_weight = 0
        for k in range(0,len(molecule.atom_list)):
                element = molecule.atom_list[k].type.split('.')[0]
                cmass[0] += molecule.atom_list[k].X * atom_mass[element]
                cmass[1] += molecule.atom_list[k].Y * atom_mass[element]
                cmass[2] += molecule.atom_list[k].Z * atom_mass[element]
                centroid[0] += molecule.atom_list[k].X
                centroid[1] += molecule.atom_list[k].Y
                centroid[2] += molecule.atom_list[k].Z
                molecular_weight += atom_mass[element]
        #print "Molecular Weight =",molecular_weight
        cmass[0] /= molecular_weight
        cmass[1] /= molecular_weight
        cmass[2] /= molecular_weight
        centroid[0] /= len(molecule.atom_list)
        centroid[1] /= len(molecule.atom_list)
        centroid[2] /= len(molecule.atom_list)
        #print 'Centroid =',centroid
        return cmass
#################################################################################################################
def molecular_weight(molecule):
        # Dictionary of atomic weights of elements
        atom_mass = {'O':15.9994 ,'N':14.00674 ,'C':12.011 ,'F':18.9984032 ,'Cl':35.4527 ,'Br':79.904
        ,'I':126.90447 ,'H':1.00794 ,'B':10.811 ,'S':32.066 ,'P':30.973762 ,'Li':6.941 ,'Na':22.98968
        ,'Mg':24.3050 ,'Al':26.981539 ,'Si':28.0855 ,'K':39.0983 ,'Ca':40.078 ,'Cr':51.9961 ,'Mn':54.93805
        ,'Fe':55.847 ,'Co':58.93320 ,'Cu':63.546 ,'Zn':65.39 ,'Se':78.96 ,'Mo':95.94 ,'Sn':118.710 ,'LP':0.0 }

        molecular_weight = 0
        for k in range(0,len(molecule.atom_list)):
                element = molecule.atom_list[k].type.split('.')[0]
                molecular_weight += atom_mass[element]
        return molecular_weight
#################################################################################################################
def calc_dipole_moment(molecule):
    uIsum=0
    uJsum=0
    uKsum=0
    dipolemoment=0
    conversion = 4.796 # Convert partialcharge*angstroms --> Coulombs*meters (Debye)

    cmass = centre_of_mass(molecule)
    #print "Centre of mass = ",cmass

    #cmass = [molecule.atom_list[0].X, molecule.atom_list[0].Y, molecule.atom_list[0].Z]
    for k in range(0,len(molecule.atom_list)):
        uIsum += molecule.atom_list[k].Q * (molecule.atom_list[k].X - cmass[0])
        uJsum += molecule.atom_list[k].Q * (molecule.atom_list[k].Y - cmass[1])
        uKsum += molecule.atom_list[k].Q * (molecule.atom_list[k].Z - cmass[2])

    umag          = sqrt( (uIsum*uIsum) + (uJsum*uJsum) + (uKsum*uKsum) )
    dipolemoment  = umag*conversion;
    uvector = [uIsum,uJsum,uKsum]

    return uvector, dipolemoment
#################################################################################################################
# Takes a single Mol object and returns a Mol object without the hydrogens
# Have to remove H from atom_list, bond_list and residue_list
def remove_hydrogens(m):
    atom_list = []
    bond_list = []
    residue_list = {}

    # Retain only heavy atoms in atom_list
    num_hvy_atoms = 0
    for i in range(len(m.atom_list)):
        if (m.atom_list[i].heavy_atom):
           atom_list.append(m.atom_list[i])
           num_hvy_atoms+=1

    # Retain only bonds containing heavy atoms
    for bond_id in range(len(m.bond_list)):
        retain_bond = True
        for atom_id in range(len(m.atom_list)):
           if (m.atom_list[atom_id].heavy_atom):
              continue  
	   # Atoms down here are always hydrogen 
           if (m.bond_list[bond_id].a1_num == m.atom_list[atom_id].num):
              retain_bond = False 
           if (m.bond_list[bond_id].a2_num == m.atom_list[atom_id].num):
              retain_bond = False
        if (retain_bond):
            bond_list.append(m.bond_list[bond_id])

    # Assuming that residue list does not change

    data = Mol(m.name,atom_list,bond_list,m.residue_list)
    ID_heavy_atoms(data.atom_list);
    return data
#################################################################################################################

