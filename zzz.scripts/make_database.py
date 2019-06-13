#! /usr/bin/python
import math, sys
import os.path
import sys
import os

#################################################################################################################
## This code is used to sort a multimol2 file by rotatable bonds and then break the resulting multimol2 into chunks
## This is used to perpare ZINC for docking.
## By William T. Berger,  
## Modified by Trent Balius
#################################################################################################################

def getrotinputs(listname,minrotsize,maxrotsize):
        os.system("pwd")
        dockpath = os.getenv ("DOCKHOME")
        dockparam = dockpath + "/parameters"

	#count_rot = 0
	count_rot = minrotsize
	while (count_rot <= maxrotsize):
                # Write dock input file 
		file = open("make_rotatable_input_"+str(count_rot)+".in",'w')
		file.write("ligand_atom_file "+listname +"\n")
		file.write("limit_max_ligands no\n")
		file.write("skip_molecule no\n")
		file.write("read_mol_solvation no\n")
		file.write("calculate_rmsd no\n")
		file.write("use_database_filter yes\n")
		file.write("dbfilter_max_heavy_atoms 9999\n")
		file.write("dbfilter_min_heavy_atoms 0\n")
		file.write("dbfilter_max_rot_bonds "+str(count_rot)+"\n")
		file.write("dbfilter_min_rot_bonds "+str(count_rot)+"\n")
		file.write("dbfilter_max_molwt 9999\n")
		file.write("dbfilter_min_molwt 0\n")
		file.write("dbfilter_max_formal_charge 5.0\n")
		file.write("dbfilter_min_formal_charge -5.0\n")
		file.write("orient_ligand no\n")
		file.write("use_internal_energy no\n")
		file.write("flexible_ligand no\n")
		file.write("bump_filter no\n")
		file.write("score_molecules no\n")
		file.write("atom_model all\n")
		file.write("vdw_defn_file "+ dockparam+"/vdw_AMBER_parm99.defn\n")
		file.write("flex_defn_file "+ dockparam+"/flex.defn\n")
		file.write("flex_drive_file "+ dockparam+"/flex_drive.tbl\n")
		file.write("ligand_outfile_prefix "+str(count_rot)+"\n")
		file.write("write_orientations no\n")
		file.write("num_scored_conformers 1\n")
		file.write("rank_ligands no\n")
                file.close()
                # run dock to brake up mol2 file by rotatable bonds.  each rotatable bond number will have a one mol2 file
		dockrun = dockpath + "/bin/dock6"
		os.system(dockrun + " -i "+ "make_rotatable_input_"+str(count_rot)+".in" + " -o " "make_rotatable_output_"+str(count_rot)+".out")
                # combine mol2 back in to one mol2 sorted by rotatble bonds.
		os.system("cat " +str(count_rot)+"_scored.mol2 >> temp.mol2")
		os.system("rm "+str(count_rot)+"_scored.mol2")
		os.system("rm make_rotatable_input_*.in")
		count_rot = count_rot + 1
	return

def getchunks(chunksize,libsize,chunknum):
        os.system("pwd")
        # chunknum is the starting lable for the file names 
        num_chunk = libsize / chunksize
	print "number of chunks: "+str(num_chunk)
	dockpath = os.getenv ("DOCKHOME")
        dockparam = dockpath + "/parameters"

	chunk_max = chunksize
	chunk_min = 0
	
	chunk = chunknum
        while (chunk <= (num_chunk + chunknum)):
		if (chunk_min == 0):
			flag = ("no")
		else: 
			flag = ("yes")
                # write dock input file
		file = open("make_chunk_input_"+str(chunk)+".in",'w')
		file.write("ligand_atom_file temp.mol2\n")
                file.write("limit_max_ligands yes\n")
                file.write("max_ligands "+str(chunk_max)+"\n")
		file.write("skip_molecule "+str(flag)+"\n")
                file.write("initial_skip "+str(chunk_min)+"\n")
		file.write("read_mol_solvation no\n")
                file.write("calculate_rmsd no\n")
                file.write("use_database_filter no\n")
                file.write("orient_ligand no\n")
                file.write("use_internal_energy no\n")
                file.write("flexible_ligand no\n")
                file.write("bump_filter no\n")
                file.write("score_molecules no\n")
                file.write("ligand_outfile_prefix chunk"+str(chunk)+"\n")
		file.write("write_orientations no\n")
                file.write("num_scored_conformers 1\n")
                file.write("rank_ligands no\n")
                file.close()
                # run dock to brake up into equal chucks.  we my onto scale chuncks (few RB) have more molecules.
                dockrun = dockpath + "/bin/dock6"
                os.system(dockrun + " -i "+ "make_chunk_input_"+str(chunk)+".in" + " -o " "make_chunk_output_"+str(chunk)+".out")
                chunk = chunk + 1
		chunk_min = chunk_min + chunksize
		chunk_max = chunk_max + chunksize
        os.system("rm temp.mol2")
	os.system("rm make_chunk_input*.in")
	return

#################################################################################################################
def num_mol_mol2(mol2file):

    fin = open(mol2file,'r')
    lines = fin.readlines()
    fin.close()

    count = 0
    for line in lines:
     linestrip = line.strip()
     if linestrip == "@<TRIPOS>MOLECULE":
        count = count + 1
   
    return count

#################################################################################################################
def main():
    if (len(sys.argv) != 7 ): # if no input
       print " This script needs the following:"
       print " (1) ZINC database multimol2 "
       print " (2) Min rot bonds "
       print " (3) Max rot bonds "
       print " (4) chunk size "
       print " (5) number to start chunk numbering "
       print " (6) library size "
       print " syntax: python make_database.py ZINCTOTAL.mol2 0 15 100000 800000" 
       return

    listname = str(sys.argv[1])
    minrotsize = int(sys.argv[2])
    maxrotsize = int(sys.argv[3])
    chunksize = int(sys.argv[4])
    chunknum = int(sys.argv[5])
    libsize = int(sys.argv[6])
    print "file input name "+listname
    print "min rot bonds "+str(minrotsize)
    print "max rot bonds "+str(maxrotsize)
    print "chunk size "+str(chunksize)
    print "chunk num "+str(chunknum)
    print "library size "+str(libsize)
    
    getrotinputs(listname,minrotsize,maxrotsize)

    num_mol_in_range =  num_mol_mol2('temp.mol2')
    print "library size in range "+str(num_mol_in_range)
    
    getchunks(chunksize,num_mol_in_range,chunknum)
    return 

#################################################################################################################
main()
