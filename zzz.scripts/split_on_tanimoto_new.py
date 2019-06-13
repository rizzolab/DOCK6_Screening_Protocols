### This scripts takes denovo mol2 with tanimoto values and separates them into temp.split.${FNUM}.mol2
### where each file contains all molecules that have the same tanimoto value to a given reference
### It also takes two inputs sys.argv[1] is the dock mol2 file this script will split up
### sys.argv[2] (if used, otherwise defaults to zero) is implemented just in case

import sys 
start_num = 0
if (len(sys.argv) > 1):
    inpt = sys.argv[1]
else:
    inpt = 'rescore_scored.mol2'
if (len(sys.argv) > 2):
    start_num = int(sys.argv[2])

with open(inpt, 'r') as fil:
    mol = ""
    tan_vals = []
    fil_dest = 0
    first_line = True
    ### reads in lines of the file until we have a whole molecule saved in mol
    for line in fil:
        if (( 'Name:' in line ) and ( not first_line  ) ):
            with open('temp.split.'+str(start_num + fil_dest)+'.mol2', 'a') as out:
		### by the time we reach this condition mol will contain a whole mol2 record 
		### and fil_dest will be set
                out.write(mol)
            mol = ""
        mol = mol + line
        if (first_line): ### this keeps the program from writing an empty file 
            first_line = False
        if ( 'Tanimoto' in line ):
            spl = line.split()
            tan = spl[2].rstrip()

            seen_val = False
            for i in range( len(tan_vals)):
                if ( tan == tan_vals[i] ): ### we've already seen a molecule with this tanimoto value
                    seen_val = True       
                    fil_dest = i           ### we know where it should go, one file per "i" 
            if (not seen_val):             ### new tanimoto value
                fil_dest = len(tan_vals)   ### its destination is in file number as many as we have so far plus one
                tan_vals.append(tan)
                    
    with open('temp.split.'+str(start_num + fil_dest)+'.mol2', 'a') as out: ### catch for last molecule in file
            out.write(mol)        
